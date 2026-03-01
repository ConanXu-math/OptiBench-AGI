"""Batch execution utilities for OptiBench workflow."""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from toolkits.paper_tools import download_paper, search_arxiv


def run_batch_pipeline(args: Any, wf: Any) -> Path:
    """Execute workflow for one or many papers and write run summary."""
    if not args.query and not args.arxiv_id:
        raise ValueError("Provide either a search query or --arxiv-id.")
    if args.top_k < 1:
        raise ValueError("--top-k must be >= 1")

    search_results: list[dict[str, object]] = []
    if args.arxiv_id:
        paper_ids = [args.arxiv_id]
    else:
        papers = search_arxiv(args.query, max_results=args.top_k, domain=args.domain)
        if not papers:
            raise RuntimeError(f"No papers found for query: {args.query}")
        search_results = [
            {
                "rank": i + 1,
                "arxiv_id": p.get("arxiv_id", ""),
                "title": p.get("title", ""),
                "summary": p.get("summary", ""),
            }
            for i, p in enumerate(papers)
        ]
        paper_ids = [p["arxiv_id"] for p in papers if p.get("arxiv_id")]
        if not paper_ids:
            raise RuntimeError("Search results contained no valid arXiv IDs.")

    summary_items: list[dict[str, object]] = []
    success = 0
    failed = 0
    logger = logging.getLogger("optibench")

    for idx, arxiv_id in enumerate(paper_ids, 1):
        paper_dir = Path(args.dataset_root) / arxiv_id.replace("/", "_")
        search_stage = "skipped" if args.arxiv_id else "ok"
        print("\n" + "=" * 60)
        print(f"[{idx}/{len(paper_ids)}] Processing arXiv {arxiv_id}")
        print("=" * 60)

        try:
            dl = download_paper(
                arxiv_id=arxiv_id,
                out_dir=args.dataset_root,
                convert_to_md=True,
                save_md_file=True,
            )
            paper_md = dl.get("md_content") or ""
            paper_name = dl.get("title") or "Untitled"
        except Exception as exc:
            failed += 1
            logger.exception("Download stage failed for arXiv %s: %s", arxiv_id, exc)
            summary_items.append(
                {
                    "arxiv_id": arxiv_id,
                    "paper_name": "",
                    "status": "failed",
                    "error": f"download failed: {exc}",
                    "output_dir": str(paper_dir),
                    "stages": {
                        "search": search_stage,
                        "download": "failed",
                        "extract": "not_run",
                        "codegen": "not_run",
                        "formalize": "skipped" if args.skip_lean else "not_run",
                        "organize": "partial",
                    },
                }
            )
            continue

        try:
            item = wf.run_from_paper(
                paper_md=paper_md,
                arxiv_id=arxiv_id,
                paper_name=paper_name,
            )
            success += 1
            stage_extract = _outline_stage_status(item.outline)
            summary_items.append(
                {
                    "arxiv_id": item.arxiv_id or arxiv_id,
                    "paper_name": item.paper_name,
                    "status": "success",
                    "output_dir": f"{args.dataset_root}/{(item.arxiv_id or arxiv_id).replace('/', '_')}",
                    "stages": {
                        "search": search_stage,
                        "download": "ok" if (paper_dir / "paper.pdf").exists() or (paper_dir / "paper.md").exists() else "failed",
                        "extract": stage_extract,
                        "codegen": "ok" if (paper_dir / "solve.py").exists() else "failed",
                        "formalize": "prove_cot_only" if args.skip_lean else ("ok" if (paper_dir / "Formal.lean").exists() else "failed"),
                        "organize": "ok" if (paper_dir / "benchmark.json").exists() and (paper_dir / "outline.json").exists() else "failed",
                    },
                }
            )
            print("Pipeline complete!")
            print(f"  Paper : {item.paper_name}")
            print(f"  arXiv : {item.arxiv_id or arxiv_id}")
            print(f"  Vars  : {item.outline.variables}")
        except Exception as exc:
            failed += 1
            logger.exception("Pipeline failed for arXiv %s: %s", arxiv_id, exc)
            stage_extract = _outline_stage_status_from_file(paper_dir / "outline.json")
            summary_items.append(
                {
                    "arxiv_id": arxiv_id,
                    "paper_name": "",
                    "status": "failed",
                    "error": str(exc),
                    "output_dir": str(paper_dir),
                    "stages": {
                        "search": search_stage,
                        "download": "ok" if (paper_dir / "paper.pdf").exists() or (paper_dir / "paper.md").exists() else "failed",
                        "extract": stage_extract,
                        "codegen": "ok" if (paper_dir / "solve.py").exists() else "not_run",
                        "formalize": "prove_cot_only" if args.skip_lean else ("ok" if (paper_dir / "Formal.lean").exists() else "not_run"),
                        "organize": "ok" if (paper_dir / "benchmark.json").exists() and (paper_dir / "outline.json").exists() else "partial",
                    },
                }
            )

    dataset_root = Path(args.dataset_root)
    dataset_root.mkdir(parents=True, exist_ok=True)
    summary = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "query": args.query,
        "requested_top_k": args.top_k,
        "total": len(paper_ids),
        "success": success,
        "failed": failed,
        "search_results": search_results,
        "items": summary_items,
    }
    summary_path = dataset_root / "run_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    if search_results:
        (dataset_root / "search_results.json").write_text(
            json.dumps(search_results, indent=2, ensure_ascii=False), encoding="utf-8"
        )

    print("\n" + "=" * 60)
    print("Batch finished!")
    print(f"  Total   : {len(paper_ids)}")
    print(f"  Success : {success}")
    print(f"  Failed  : {failed}")
    print(f"  Summary : {summary_path}")
    print("=" * 60)
    return summary_path


def _outline_stage_status(outline: object) -> str:
    objective = getattr(outline, "objective", "") or ""
    constraints = getattr(outline, "constraints", []) or []
    variables = getattr(outline, "variables", []) or []
    if objective.strip() or constraints or variables:
        return "ok"
    return "failed"


def _outline_stage_status_from_file(path: Path) -> str:
    if not path.exists():
        return "not_run"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return "failed"
    objective = (data.get("objective") or "").strip()
    constraints = data.get("constraints") or []
    variables = data.get("variables") or []
    if objective or constraints or variables:
        return "ok"
    return "failed"

