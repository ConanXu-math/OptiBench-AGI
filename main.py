"""OptiBench-AGI – CLI entry point.

Usage:
    # Basic search query
    python main.py "convex optimization relaxation techniques"

    # Specific arXiv paper ID
    python main.py --arxiv-id 2301.12345

    # Custom output directory & retries
    python main.py "sparse recovery" --dataset-root ./my_dataset --max-retries 5
"""

from __future__ import annotations

import argparse
import logging
import os
import sys


def _build_model():
    """Resolve the LLM model from env vars.

    Environment variables:
        OPTIBENCH_PROVIDER   openai | openai-like | anthropic | google
        OPTIBENCH_MODEL      model id  (default: gpt-4o)
        OPTIBENCH_BASE_URL   custom endpoint for openai-like provider
        OPTIBENCH_API_KEY    API key (falls back to provider-specific vars)
    """
    provider = os.getenv("OPTIBENCH_PROVIDER", "openai").lower()
    model_id = os.getenv("OPTIBENCH_MODEL", "gpt-4o")
    base_url = os.getenv("OPTIBENCH_BASE_URL", "")
    api_key = os.getenv("OPTIBENCH_API_KEY", "")

    if provider in ("openai-like", "openailike"):
        from agno.models.openai.like import OpenAILike

        kwargs: dict = {"id": model_id}
        if base_url:
            kwargs["base_url"] = base_url
        if api_key:
            kwargs["api_key"] = api_key
        return OpenAILike(**kwargs)

    if provider == "openai":
        from agno.models.openai import OpenAIChat

        return OpenAIChat(id=model_id)

    if provider == "anthropic":
        from agno.models.anthropic import Claude

        return Claude(id=model_id)

    if provider == "google":
        from agno.models.google import Gemini

        return Gemini(id=model_id)

    raise ValueError(
        f"Unknown OPTIBENCH_PROVIDER={provider!r}. "
        "Supported: openai, openai-like, anthropic, google"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="OptiBench-AGI: automated optimisation paper benchmarking",
    )
    parser.add_argument(
        "query",
        nargs="?",
        default=None,
        help="Free-text search query for arXiv (e.g. 'convex optimisation').",
    )
    parser.add_argument(
        "--arxiv-id",
        default=None,
        help="Directly specify an arXiv paper ID instead of searching.",
    )
    parser.add_argument(
        "--dataset-root",
        default="./dataset",
        help="Root directory for saving benchmark artefacts (default: ./dataset).",
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=3,
        help="Max validation retry attempts per stage (default: 3).",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable debug-level logging."
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(name)-12s  %(levelname)-8s  %(message)s",
    )

    if not args.query and not args.arxiv_id:
        parser.error("Provide either a search query or --arxiv-id.")

    query = args.query or f"arXiv paper {args.arxiv_id}"
    if args.arxiv_id:
        query = (
            f"Download and convert the paper with arXiv ID {args.arxiv_id}. "
            f"Use download_and_convert_paper directly."
        )

    from workflow import OptiBenchWorkflow

    model = _build_model()
    wf = OptiBenchWorkflow(
        model=model,
        dataset_root=args.dataset_root,
        max_retries=args.max_retries,
    )

    item = wf.run(query)

    print("\n" + "=" * 60)
    print("Pipeline complete!")
    print(f"  Paper : {item.paper_name}")
    print(f"  arXiv : {item.arxiv_id}")
    print(f"  Vars  : {item.outline.variables}")
    print(f"  Output: {args.dataset_root}/{item.arxiv_id or 'unknown'}/")
    print("=" * 60)


if __name__ == "__main__":
    main()
