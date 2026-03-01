"""OptiBenchWorkflow – orchestrates the full benchmarking pipeline.

Pipeline stages:
  1. Search & Parse   → Research Agent
  2. Extract Outline  → Extraction Agent
  3. Code Generation  → Coding Agent  + validate_python_code (retry ≤ N)
  4. Formalisation     → Formal Agent  + validate_lean_code   (retry ≤ N)
  5. Persist           → save BenchmarkItem to disk
"""

from __future__ import annotations

import json
import logging
import os
import re
from pathlib import Path
from typing import TYPE_CHECKING

from agents.coding import create_coding_agent
from agents.extraction_critic import create_extraction_critic_agent
from agents.extraction import create_extraction_agent
from agents.formalization import create_formalization_agent
from agents.locator import create_locator_agent
from agents.research import create_research_agent
from quality.gate import evaluate_quality_gate
from schema.models import (
    BenchmarkItem,
    CodingOutput,
    ExtractionCritique,
    FormalizationOutput,
    MathOutline,
    NotationItem,
    SectionSelection,
)
from toolkits.validators import validate_lean_code, validate_python_code

if TYPE_CHECKING:
    from agno.agent import Agent
    from agno.models.base import Model

logger = logging.getLogger("optibench")


class OptiBenchWorkflow:
    """End-to-end pipeline from arXiv query to validated BenchmarkItem."""

    def __init__(
        self,
        model: Model,
        *,
        dataset_root: str | Path = "./dataset",
        max_retries: int = 3,
        skip_lean: bool = False,
        require_outline: bool = True,
    ) -> None:
        self.model = model
        self.dataset_root = Path(dataset_root)
        self.max_retries = max_retries
        self.skip_lean = skip_lean
        self.require_outline = require_outline

        self.research_agent: Agent = create_research_agent(model)
        self.locator_agent: Agent = create_locator_agent(model)
        self.extraction_agent: Agent = create_extraction_agent(model)
        self.extraction_critic_agent: Agent = create_extraction_critic_agent(model)
        self.coding_agent: Agent = create_coding_agent(model)
        self.formalization_agent: Agent = create_formalization_agent(model)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self, query: str) -> BenchmarkItem:
        """Execute the full pipeline for *query* and return a BenchmarkItem."""
        logger.info("=== Step 1: Search & Parse ===")
        paper_md, arxiv_id, paper_name = self._step_search(query)
        return self.run_from_paper(
            paper_md=paper_md,
            arxiv_id=arxiv_id,
            paper_name=paper_name,
        )

    def run_from_paper(
        self,
        *,
        paper_md: str,
        arxiv_id: str,
        paper_name: str,
    ) -> BenchmarkItem:
        """Run pipeline from an already-downloaded paper markdown.

        This method intentionally skips Step 1 so search/download can be
        orchestrated independently by the caller.
        """

        logger.info("=== Step 2: Extract MathOutline ===")
        outline = self._step_extract(paper_md)
        if not _is_outline_meaningful(outline):
            logger.warning(
                "  Outline is empty or minimal (no objective/variables). "
                "Coding step may produce generic code. Check PDF conversion and paper content."
            )
            if self.require_outline:
                raise RuntimeError(
                    "Extraction produced empty outline (use --allow-empty-outline to continue anyway). "
                    "Ensure PDF→Markdown works (marker/nougat) or try another paper."
                )

        logger.info("=== Step 3: Generate & Validate Python Code ===")
        coding_out = self._step_code(outline)

        if self.skip_lean:
            logger.info("=== Step 4: Lean 4 (skipped), generating LaTeX prove_cot ===")
            prove_cot, lean4_formal = self._step_prove_cot(outline), ""
        else:
            logger.info("=== Step 4: Generate & Validate Lean 4 Code ===")
            formal_out = self._step_formalize(outline)
            prove_cot, lean4_formal = formal_out.prove_cot, formal_out.lean4_formal

        logger.info("=== Step 5: Persist ===")
        item = BenchmarkItem(
            paper_name=paper_name,
            arxiv_id=arxiv_id,
            outline=outline,
            prove_cot=prove_cot,
            pseudocode=coding_out.pseudocode,
            pycode=coding_out.pycode,
            lean4_formal=lean4_formal,
        )
        gate = evaluate_quality_gate(item, skip_lean=self.skip_lean)
        if not gate.passed:
            logger.warning("  Quality gate failed: %s", "; ".join(gate.notes))
        else:
            logger.info("  Quality gate passed.")
        self._save(item, paper_md=paper_md, arxiv_id=arxiv_id)
        return item

    # ------------------------------------------------------------------
    # Step implementations
    # ------------------------------------------------------------------

    def _step_search(self, query: str) -> tuple[str, str, str]:
        """Step 1: search arXiv, download & convert the best paper."""
        response = self.research_agent.run(
            f"Search arXiv for the following optimisation topic and return the "
            f"full Markdown of the most relevant paper:\n\n{query}"
        )
        paper_md: str = response.content or ""

        # OpenAI-like providers may intermittently fail tool-calling for this step.
        # If that happens, fall back to deterministic local tooling.
        if _is_agent_error_text(paper_md) or len(paper_md.strip()) < 200:
            logger.warning(
                "  Research agent output is invalid/too short; "
                "falling back to direct arXiv search + download pipeline."
            )
            fallback = self._step_search_fallback(query)
            if fallback is not None:
                return fallback

        if "ERROR:" in paper_md or "ERROR " in paper_md:
            logger.warning(
                "  PDF conversion may have failed (error in content). "
                "Install marker or nougat for full text: pip install marker-pdf | nougat-ocr"
            )
        if len(paper_md) < 500:
            logger.warning("  Paper content very short (%d chars); extraction may be poor.", len(paper_md))
        arxiv_id = _extract_arxiv_id(paper_md)
        paper_name = _extract_title(paper_md)
        logger.info("  Paper: %s  (arXiv: %s)", paper_name, arxiv_id)
        return paper_md, arxiv_id, paper_name

    def _step_search_fallback(self, query: str) -> tuple[str, str, str] | None:
        """Deterministic fallback: search arXiv and download paper without LLM."""
        try:
            from toolkits.paper_tools import download_paper, search_arxiv

            results = search_arxiv(query, max_results=5)
            if not results:
                logger.error("  Fallback search found no arXiv papers for query: %s", query)
                return None

            best = results[0]
            arxiv_id = best.get("arxiv_id", "")
            if not arxiv_id:
                logger.error("  Fallback search returned paper without arXiv ID.")
                return None

            dl = download_paper(
                arxiv_id=arxiv_id,
                out_dir=self.dataset_root,
                convert_to_md=True,
                save_md_file=True,
            )
            paper_md = dl.get("md_content") or ""
            paper_name = dl.get("title") or best.get("title") or "Untitled"

            if not paper_md.strip():
                logger.error("  Fallback download succeeded but markdown content is empty.")
                return None

            logger.info("  [fallback] Paper: %s  (arXiv: %s)", paper_name, arxiv_id)
            return paper_md, arxiv_id, paper_name
        except Exception as exc:
            logger.exception("  Fallback search pipeline failed: %s", exc)
            return None

    def _step_extract(self, paper_md: str) -> MathOutline:
        """Step 2: extract a structured MathOutline from paper Markdown."""
        max_attempts = max(3, self.max_retries)
        logger.info("  Step 2a: Locator agent selecting key sections")
        focused_md = self._step_locate_sections(paper_md)
        prompt_prefix = (
            "Extract the primary optimisation formulation from the following "
            "paper Markdown.\n\n"
        )
        last_error = ""
        last_outline = MathOutline(objective="", constraints=[], variables=[], notation_table=[])
        budgets = _extraction_budgets(max_attempts)

        for attempt in range(1, max_attempts + 1):
            if attempt > 1:
                prompt_prefix = (
                    "The previous extraction failed or was empty. "
                    "You MUST return valid JSON matching MathOutline. "
                    "Return non-empty 'objective' (LaTeX) and at least 'variables' or 'constraints'. "
                    "Strictly follow the critic fix instructions if provided.\n\n"
                    "Paper Markdown:\n\n"
                )

            budget = budgets[min(attempt - 1, len(budgets) - 1)]
            extraction_md = _prepare_extraction_input(focused_md, char_budget=budget)
            logger.info(
                "  Extraction input length: %d chars (attempt %d/%d, budget=%d)",
                len(extraction_md),
                attempt,
                max_attempts,
                budget,
            )
            response = self.extraction_agent.run(f"{prompt_prefix}{extraction_md}")
            outline, error = _parse_outline_output(response.content)
            if outline is None:
                last_error = error or "unknown extraction error"
                logger.warning(
                    "  Extraction attempt %d/%d failed: %s",
                    attempt,
                    max_attempts,
                    last_error[:200],
                )
                continue

            last_outline = outline
            logger.info(
                "  Objective: %s  |  %d constraints  |  %d variables",
                (outline.objective[:80] + "…") if len(outline.objective) > 80 else outline.objective or "(empty)",
                len(outline.constraints),
                len(outline.variables),
            )
            critique = self._step_critique_extraction(extraction_md, outline)
            if _is_outline_meaningful(outline) and critique.is_acceptable:
                return outline
            if critique.fix_prompt.strip():
                prompt_prefix = (
                    "The previous extraction was rejected by quality critic.\n"
                    f"Fix instruction:\n{critique.fix_prompt.strip()}\n\n"
                    "Paper Markdown:\n\n"
                )
            logger.warning(
                "  Outline rejected (acceptable=%s, score=%d), retrying extraction.",
                critique.is_acceptable,
                critique.score,
            )

        if _is_outline_meaningful(last_outline):
            return last_outline
        if _is_timeout_like_error(last_error):
            logger.warning(
                "  Extraction timed out repeatedly. Using deterministic fallback outline."
            )
            return _fallback_outline_from_markdown(paper_md)
        raise RuntimeError(
            "Extraction failed after retries. "
            f"Last error: {last_error or 'empty outline returned repeatedly'}"
        )

    def _step_locate_sections(self, paper_md: str) -> str:
        """Use locator agent to focus extraction context."""
        if len(paper_md) < 5000:
            return paper_md
        prompt = (
            "Select high-value sections for extracting optimization formulation.\n\n"
            f"{paper_md}"
        )
        try:
            response = self.locator_agent.run(prompt)
            raw = response.content
        except Exception as e:
            logger.warning(
                "  Locator API failed (%s); using heuristic extraction input.",
                type(e).__name__,
            )
            return _prepare_extraction_input(paper_md, char_budget=12000)
        if isinstance(raw, SectionSelection):
            selected = raw.selected_markdown.strip()
            if selected:
                logger.info("  Locator selected %d chars.", len(selected))
                return selected
        elif isinstance(raw, str) and raw.strip():
            try:
                parsed = SectionSelection.model_validate_json(raw)
                if parsed.selected_markdown.strip():
                    logger.info("  Locator selected %d chars.", len(parsed.selected_markdown))
                    return parsed.selected_markdown
            except Exception:
                pass
        logger.warning("  Locator failed; using heuristic extraction input.")
        return _prepare_extraction_input(paper_md, char_budget=12000)

    def _step_critique_extraction(
        self,
        extraction_md: str,
        outline: MathOutline,
    ) -> ExtractionCritique:
        """Critique extraction quality and produce fix prompt."""
        prompt = (
            "Review extraction quality.\n\n"
            f"Source excerpt:\n```markdown\n{extraction_md[:6000]}\n```\n\n"
            f"Extracted outline:\n```json\n{outline.model_dump_json(indent=2)}\n```"
        )
        response = self.extraction_critic_agent.run(prompt)
        raw = response.content
        if isinstance(raw, ExtractionCritique):
            return raw
        if isinstance(raw, str):
            try:
                if "```json" in raw:
                    raw = raw.split("```json", 1)[1].split("```", 1)[0].strip()
                return ExtractionCritique.model_validate_json(raw)
            except Exception:
                pass
        # conservative default: request retry when critic unavailable
        return ExtractionCritique(
            score=30,
            is_acceptable=False,
            issues=["Critic unavailable or malformed output."],
            fix_prompt="Ensure objective, constraints and variables strictly match source equations.",
        )

    def _step_code(self, outline: MathOutline) -> CodingOutput:
        """Step 3: generate Python code with validation retry loop."""
        outline_json = outline.model_dump_json(indent=2)
        prompt = (
            f"Generate Python code and pseudocode for the following optimisation "
            f"problem:\n\n```json\n{outline_json}\n```"
        )

        last_output: CodingOutput | None = None
        for attempt in range(1, self.max_retries + 1):
            logger.info("  Coding attempt %d/%d", attempt, self.max_retries)
            response = self.coding_agent.run(prompt)
            last_output = _parse_coding_output(response.content)
            if last_output is None:
                logger.warning("  Coding agent returned invalid output, retrying.")
                prompt = (
                    f"The previous response was invalid (API error or malformed output). "
                    f"Please try again and return valid JSON with 'pseudocode' and 'pycode'.\n\n"
                    f"Original MathOutline:\n```json\n{outline_json}\n```"
                )
                continue

            validation = validate_python_code(last_output.pycode)
            if validation.startswith("SUCCESS"):
                logger.info("  Python validation passed.")
                return last_output

            logger.warning("  Validation failed:\n%s", validation)
            prompt = (
                f"The previous Python code FAILED validation.\n\n"
                f"Error:\n```\n{validation}\n```\n\n"
                f"Previous code:\n```python\n{last_output.pycode}\n```\n\n"
                f"Original MathOutline:\n```json\n{outline_json}\n```\n\n"
                f"Fix the code so it runs without errors."
            )

        logger.error("  Exhausted %d coding attempts.", self.max_retries)
        if last_output is None:
            raise RuntimeError(
                "Coding agent did not return valid output (API errors or malformed response). "
                "Check your API key and try again."
            )
        return last_output

    def _step_formalize(self, outline: MathOutline) -> FormalizationOutput:
        """Step 4: generate Lean 4 code with validation retry loop."""
        outline_json = outline.model_dump_json(indent=2)
        prompt = (
            f"Formalize the following optimisation problem in Lean 4 with Mathlib:\n\n"
            f"```json\n{outline_json}\n```"
        )

        last_output: FormalizationOutput | None = None
        for attempt in range(1, self.max_retries + 1):
            logger.info("  Lean 4 attempt %d/%d", attempt, self.max_retries)
            response = self.formalization_agent.run(prompt)
            last_output = _parse_formalization_output(response.content)
            if last_output is None:
                logger.warning("  Formalization agent returned invalid output, retrying.")
                prompt = (
                    f"The previous response was invalid (API error or malformed output). "
                    f"Please try again and return valid JSON with 'prove_cot' and 'lean4_formal'.\n\n"
                    f"Original MathOutline:\n```json\n{outline_json}\n```"
                )
                continue

            validation = validate_lean_code(last_output.lean4_formal)
            if validation.startswith("SUCCESS"):
                logger.info("  Lean 4 validation passed.")
                return last_output

            logger.warning("  Lean validation failed:\n%s", validation)
            prompt = (
                f"The previous Lean 4 code FAILED compilation.\n\n"
                f"Error:\n```\n{validation}\n```\n\n"
                f"Previous code:\n```lean\n{last_output.lean4_formal}\n```\n\n"
                f"Original MathOutline:\n```json\n{outline_json}\n```\n\n"
                f"Fix the Lean 4 code so it compiles with `lake build`."
            )

        logger.error("  Exhausted %d formalisation attempts.", self.max_retries)
        if last_output is None:
            raise RuntimeError(
                "Formalization agent did not return valid output (API errors or malformed response)."
            )
        return last_output

    def _step_prove_cot(self, outline: MathOutline) -> str:
        """Generate LaTeX proof-style derivation without Lean formalization."""
        outline_json = outline.model_dump_json(indent=2)
        prompt = (
            "Produce ONLY a concise LaTeX derivation (prove_cot) for the optimization problem. "
            "Do NOT include Lean code. Set `lean4_formal` to an empty string.\n\n"
            f"MathOutline:\n```json\n{outline_json}\n```"
        )
        last_err = ""
        for attempt in range(1, self.max_retries + 1):
            logger.info("  prove_cot attempt %d/%d", attempt, self.max_retries)
            response = self.formalization_agent.run(prompt)
            parsed = _parse_formalization_output(response.content)
            if parsed is not None and parsed.prove_cot.strip():
                return parsed.prove_cot.strip()
            last_err = str(response.content)[:240]
            prompt = (
                "Previous response was invalid. Return valid JSON for FormalizationOutput "
                "with non-empty `prove_cot` in LaTeX and empty `lean4_formal`.\n\n"
                f"MathOutline:\n```json\n{outline_json}\n```"
            )
        logger.warning("  prove_cot generation failed, using deterministic fallback. Last: %s", last_err)
        return _fallback_prove_cot(outline)

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def _save(
        self,
        item: BenchmarkItem,
        *,
        paper_md: str | None = None,
        arxiv_id: str = "",
    ) -> Path:
        """Write all artefacts under ``dataset_root/<paper_id>/``."""
        folder_name = item.arxiv_id.replace("/", "_") if item.arxiv_id else "unknown"
        out_dir = self.dataset_root / folder_name
        out_dir.mkdir(parents=True, exist_ok=True)

        if paper_md and paper_md.strip():
            (out_dir / "paper.md").write_text(paper_md, encoding="utf-8")
            logger.info("  Saved paper.md (full text from search)")

        if arxiv_id and not (out_dir / "paper.pdf").exists():
            try:
                from toolkits.paper_tools import download_arxiv_pdf
                pdf_path = download_arxiv_pdf(arxiv_id, out_dir)
                if pdf_path.exists():
                    import shutil
                    dest = out_dir / "paper.pdf"
                    if pdf_path != dest:
                        shutil.move(str(pdf_path), str(dest))
                    logger.info("  Saved paper.pdf")
            except Exception as e:
                logger.warning("  Could not save PDF: %s", e)

        (out_dir / "benchmark.json").write_text(
            item.model_dump_json(indent=2), encoding="utf-8"
        )
        (out_dir / "outline.json").write_text(
            item.outline.model_dump_json(indent=2), encoding="utf-8"
        )
        if item.pycode:
            (out_dir / "solve.py").write_text(item.pycode, encoding="utf-8")
        if item.pseudocode:
            (out_dir / "pseudocode.txt").write_text(item.pseudocode, encoding="utf-8")
        if item.lean4_formal:
            (out_dir / "Formal.lean").write_text(item.lean4_formal, encoding="utf-8")

        logger.info("  Artefacts saved to %s", out_dir)
        return out_dir


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------


def _is_outline_meaningful(outline: MathOutline) -> bool:
    """True if outline has at least an objective or variables/constraints."""
    if outline.objective and outline.objective.strip():
        return True
    if outline.variables or outline.constraints:
        return True
    if outline.notation_table:
        return True
    return False


def _parse_coding_output(content: object) -> CodingOutput | None:
    """Parse agent response into CodingOutput; return None if invalid or API error."""
    if isinstance(content, CodingOutput):
        return content
    if not isinstance(content, str):
        return None
    raw = content.strip()
    if not raw or "Error from" in raw or "not set" in raw or "Error in Agent" in raw:
        return None
    try:
        if "```json" in raw:
            raw = raw.split("```json", 1)[1].split("```", 1)[0].strip()
        elif "```" in raw:
            raw = raw.split("```", 1)[1].split("```", 1)[0].strip()
        return CodingOutput.model_validate_json(raw)
    except Exception:
        return None


def _parse_outline_output(content: object) -> tuple[MathOutline | None, str]:
    """Parse extraction-agent response into MathOutline.

    Returns (outline, error_message). On success, error_message is empty.
    """
    if isinstance(content, MathOutline):
        return content, ""
    if not isinstance(content, str):
        return None, f"unexpected content type: {type(content).__name__}"

    raw = content.strip()
    if not raw:
        return None, "empty response"

    lowered = raw.lower()
    transient_error_patterns = [
        "request timed out",
        "api connection error",
        "error in agent run",
        "error from openai api",
        "rate limit",
        "service unavailable",
    ]
    if any(p in lowered for p in transient_error_patterns):
        return None, raw
    if raw.startswith("ERROR") or "not set" in lowered or ("api" in raw and "key" in lowered):
        return None, raw

    try:
        if "```json" in raw:
            raw = raw.split("```json", 1)[1].split("```", 1)[0].strip()
        elif "```" in raw:
            raw = raw.split("```", 1)[1].split("```", 1)[0].strip()
        return MathOutline.model_validate_json(raw), ""
    except Exception as exc:
        return None, f"invalid json: {exc}"


def _parse_formalization_output(content: object) -> FormalizationOutput | None:
    """Parse agent response into FormalizationOutput; return None if invalid or API error."""
    if isinstance(content, FormalizationOutput):
        return content
    if not isinstance(content, str):
        return None
    raw = content.strip()
    if not raw or "Error from" in raw or "not set" in raw or "Error in Agent" in raw:
        return None
    try:
        if "```json" in raw:
            raw = raw.split("```json", 1)[1].split("```", 1)[0].strip()
        elif "```" in raw:
            raw = raw.split("```", 1)[1].split("```", 1)[0].strip()
        return FormalizationOutput.model_validate_json(raw)
    except Exception:
        return None


def _is_agent_error_text(text: str) -> bool:
    if not text:
        return True
    lowered = text.lower()
    patterns = [
        "error from openai api",
        "error in agent run",
        "api key",
        "request id:",
        "error from",
        "tool call",
    ]
    return any(p in lowered for p in patterns)


_ARXIV_ID_RE = re.compile(r"(\d{4}\.\d{4,5}(?:v\d+)?)")


def _extract_arxiv_id(text: str) -> str:
    m = _ARXIV_ID_RE.search(text)
    return m.group(1) if m else ""


def _extract_title(md: str) -> str:
    for line in md.splitlines():
        stripped = line.strip()
        if stripped.startswith("# "):
            return stripped[2:].strip()
    return "Untitled"


def _extraction_budgets(max_attempts: int) -> list[int]:
    """Build decreasing char budgets for extraction retries."""
    base = int(os.getenv("OPTIBENCH_EXTRACT_MAX_CHARS", "24000"))
    levels = [base, max(base // 2, 12000), max(base // 3, 8000), 6000]
    return levels[: max(1, max_attempts)]


def _prepare_extraction_input(paper_md: str, *, char_budget: int) -> str:
    """Compress long paper markdown before sending to LLM.

    Strategy:
    1) Keep title/abstract and lines likely containing formulations.
    2) If still too short, append head/tail snippets as backup context.
    3) Cap total chars by `char_budget`.
    """
    if len(paper_md) <= char_budget:
        return paper_md

    lines = paper_md.splitlines()
    kept: list[str] = []
    key_patterns = [
        "abstract",
        "problem",
        "formulation",
        "objective",
        "constraint",
        "min ",
        "max ",
        "\\min",
        "\\max",
        "s.t.",
        "subject to",
        "lagrangian",
        "kkt",
    ]

    # Always keep early metadata block.
    kept.extend(lines[:120])

    for ln in lines:
        low = ln.lower()
        if any(p in low for p in key_patterns):
            kept.append(ln)
        if len("\n".join(kept)) >= char_budget:
            break

    merged = "\n".join(kept)
    if len(merged) < min(6000, char_budget):
        head = paper_md[: char_budget // 2]
        tail = paper_md[-(char_budget // 2) :]
        merged = head + "\n\n<!-- tail context -->\n\n" + tail

    return merged[:char_budget]


def _is_timeout_like_error(error_text: str) -> bool:
    low = (error_text or "").lower()
    return (
        "request timed out" in low
        or "api connection error" in low
        or "connection timed out" in low
        or "read timed out" in low
    )


def _fallback_outline_from_markdown(paper_md: str) -> MathOutline:
    """Heuristic extraction when LLM calls keep timing out."""
    low = paper_md.lower()

    # Try to infer algorithm-specific variable naming.
    variables = ["x"]
    if "admm" in low:
        variables = ["x", "z", "u"]
    elif "primal-dual" in low or "pdhg" in low or "chambolle" in low:
        variables = ["x", "y"]
    elif "proximal point" in low or "ppa" in low:
        variables = ["x_k"]

    has_constraints = ("subject to" in low) or ("s.t." in low) or ("constraint" in low)
    constraints = [r"c_i(x) \le 0,\; i=1,\dots,m"] if has_constraints else []
    notation = [
        NotationItem(symbol=r"x", dimension=r"\mathbb{R}^n", description="decision variable"),
        NotationItem(symbol=r"f(x)", dimension=r"\mathbb{R}", description="objective function"),
    ]
    if "admm" in low:
        objective = r"\min_{x,z}\; f(x)+g(z)\;\text{s.t.}\;Ax+Bz=c"
        notation.extend(
            [
                NotationItem(symbol=r"z", dimension=r"\mathbb{R}^m", description="auxiliary variable"),
                NotationItem(symbol=r"u", dimension=r"\mathbb{R}^p", description="scaled dual variable"),
            ]
        )
        constraints = [r"Ax + Bz = c"]
    elif "primal-dual" in low or "pdhg" in low:
        objective = r"\min_x \max_y\; \langle Kx,y\rangle + G(x)-F^*(y)"
    else:
        objective = r"\min_x f(x)"

    return MathOutline(
        objective=objective,
        constraints=constraints,
        variables=variables,
        notation_table=notation,
    )


def _fallback_prove_cot(outline: MathOutline) -> str:
    """Deterministic LaTeX derivation used when LLM prove_cot fails."""
    objective = outline.objective.strip() or r"\min_x f(x)"
    if outline.constraints:
        cons = r",\; ".join(c.strip() for c in outline.constraints if c.strip())
    else:
        cons = r"\text{(no explicit constraints)}"
    return (
        r"\textbf{Problem.}\quad "
        + objective
        + r"\;\text{s.t.}\;"
        + cons
        + "\n\n"
        + r"\textbf{Lagrangian.}\quad \mathcal{L}(x,\lambda)=f(x)+\lambda^\top g(x)."
        + "\n\n"
        + r"\textbf{First-order condition.}\quad \nabla_x \mathcal{L}(x^\star,\lambda^\star)=0,"
        + r"\; g(x^\star)\preceq 0,\; \lambda^\star \succeq 0,\; "
        + r"\lambda_i^\star g_i(x^\star)=0."
        + "\n\n"
        + r"\textbf{Conclusion.}\quad "
        + r"(x^\star,\lambda^\star)\ \text{satisfies KKT conditions and yields a stationary optimum under standard regularity assumptions.}"
    )
