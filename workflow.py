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
import re
from pathlib import Path
from typing import TYPE_CHECKING

from agents.coding import create_coding_agent
from agents.extraction import create_extraction_agent
from agents.formalization import create_formalization_agent
from agents.research import create_research_agent
from schema.models import BenchmarkItem, CodingOutput, FormalizationOutput, MathOutline
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
    ) -> None:
        self.model = model
        self.dataset_root = Path(dataset_root)
        self.max_retries = max_retries

        self.research_agent: Agent = create_research_agent(model)
        self.extraction_agent: Agent = create_extraction_agent(model)
        self.coding_agent: Agent = create_coding_agent(model)
        self.formalization_agent: Agent = create_formalization_agent(model)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self, query: str) -> BenchmarkItem:
        """Execute the full pipeline for *query* and return a BenchmarkItem."""
        logger.info("=== Step 1: Search & Parse ===")
        paper_md, arxiv_id, paper_name = self._step_search(query)

        logger.info("=== Step 2: Extract MathOutline ===")
        outline = self._step_extract(paper_md)

        logger.info("=== Step 3: Generate & Validate Python Code ===")
        coding_out = self._step_code(outline)

        logger.info("=== Step 4: Generate & Validate Lean 4 Code ===")
        formal_out = self._step_formalize(outline)

        logger.info("=== Step 5: Persist ===")
        item = BenchmarkItem(
            paper_name=paper_name,
            arxiv_id=arxiv_id,
            outline=outline,
            prove_cot=formal_out.prove_cot,
            pseudocode=coding_out.pseudocode,
            pycode=coding_out.pycode,
            lean4_formal=formal_out.lean4_formal,
        )
        self._save(item)
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

        arxiv_id = _extract_arxiv_id(paper_md)
        paper_name = _extract_title(paper_md)
        logger.info("  Paper: %s  (arXiv: %s)", paper_name, arxiv_id)
        return paper_md, arxiv_id, paper_name

    def _step_extract(self, paper_md: str) -> MathOutline:
        """Step 2: extract a structured MathOutline from paper Markdown."""
        response = self.extraction_agent.run(
            f"Extract the primary optimisation formulation from the following "
            f"paper Markdown:\n\n{paper_md}"
        )
        outline: MathOutline = response.content
        logger.info(
            "  Objective: %s  |  %d constraints  |  %d variables",
            outline.objective[:80],
            len(outline.constraints),
            len(outline.variables),
        )
        return outline

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
            last_output = response.content

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
        assert last_output is not None
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
            last_output = response.content

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
        assert last_output is not None
        return last_output

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def _save(self, item: BenchmarkItem) -> Path:
        """Write all artefacts under ``dataset_root/<paper_id>/``."""
        folder_name = item.arxiv_id.replace("/", "_") if item.arxiv_id else "unknown"
        out_dir = self.dataset_root / folder_name
        out_dir.mkdir(parents=True, exist_ok=True)

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
