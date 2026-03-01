"""OptiBench-AGI CLI entrypoint.

Architecture:
- `agents/` : LLM personas and instructions
- `toolkits/` : deterministic tools (search/download/validate/convert)
- `workflow.py` : pipeline orchestration logic
- `runner/` : batch execution and summary writing
- `cli/` : interactive prompts and utility subcommands
"""

from __future__ import annotations

import argparse
import logging
import os
import sys

from cli.commands import (
    cmd_convert_pdf,
    cmd_download,
    cmd_info,
    cmd_list,
    cmd_search,
)
from cli.interactive import interactive_pipeline_args
from runner.pipeline_batch import run_batch_pipeline

# Load .env so OPENAI_API_KEY / OPTIBENCH_* are available
try:
    from dotenv import load_dotenv

    load_dotenv()
except ImportError:
    pass


def _build_model():
    """Resolve model from env vars in one place."""
    provider = os.getenv("OPTIBENCH_PROVIDER", "openai").lower()
    model_id = os.getenv("OPTIBENCH_MODEL", "gpt-4o")
    base_url = os.getenv("OPTIBENCH_BASE_URL", "")
    api_key = os.getenv("OPTIBENCH_API_KEY", "")
    timeout_s = os.getenv("OPTIBENCH_API_TIMEOUT", "").strip()

    if provider in ("openai-like", "openailike"):
        from agno.models.openai.like import OpenAILike

        kwargs: dict = {"id": model_id}
        if base_url:
            kwargs["base_url"] = base_url
        if api_key:
            kwargs["api_key"] = api_key
        if timeout_s:
            try:
                kwargs["timeout"] = float(timeout_s)
            except ValueError:
                pass
        try:
            return OpenAILike(**kwargs)
        except TypeError:
            kwargs.pop("timeout", None)
            return OpenAILike(**kwargs)

    if provider == "openai":
        from agno.models.openai import OpenAIChat

        kwargs: dict = {"id": model_id}
        if timeout_s:
            try:
                kwargs["timeout"] = float(timeout_s)
            except ValueError:
                pass
        try:
            return OpenAIChat(**kwargs)
        except TypeError:
            kwargs.pop("timeout", None)
            return OpenAIChat(**kwargs)

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


def _dispatch_tool_subcommand() -> bool:
    """Handle utility subcommands; return True if handled."""
    if len(sys.argv) < 2:
        return False
    sub = sys.argv[1]
    if sub not in ("download", "search", "convert-pdf", "list", "info"):
        return False

    subparser = argparse.ArgumentParser(prog=f"main.py {sub}")
    if sub == "search":
        subparser.add_argument("query", help="Search keywords")
        subparser.add_argument("--max", type=int, default=10, help="Max results (default 10)")
        subparser.add_argument(
            "--domain",
            default=os.getenv("OPTIBENCH_DOMAIN", "continuous"),
            choices=["continuous", "all"],
            help="Search domain filter (default: continuous).",
        )
        args = subparser.parse_args(sys.argv[2:])
        cmd_search(args.query, args.max, args.domain)
    elif sub == "download":
        subparser.add_argument("arxiv_id", help="arXiv ID (e.g. 1406.0899 or 1406.0899v4)")
        subparser.add_argument("--output", "-o", default="./dataset", help="Output directory (default ./dataset)")
        subparser.add_argument("--no-convert", action="store_true", help="Only download PDF, do not convert to Markdown")
        args = subparser.parse_args(sys.argv[2:])
        cmd_download(args.arxiv_id, args.output, args.no_convert)
    elif sub == "list":
        subparser.add_argument("--dataset-root", "-d", default="./dataset", help="Dataset directory (default ./dataset)")
        args = subparser.parse_args(sys.argv[2:])
        cmd_list(args.dataset_root)
    elif sub == "info":
        subparser.add_argument("arxiv_id", help="arXiv ID (e.g. 1406.0899)")
        args = subparser.parse_args(sys.argv[2:])
        cmd_info(args.arxiv_id)
    else:  # convert-pdf
        subparser.add_argument("pdf_path", help="Path to PDF file")
        subparser.add_argument("-o", "--output", default=None, help="Output .md file (default: print to stdout)")
        args = subparser.parse_args(sys.argv[2:])
        cmd_convert_pdf(args.pdf_path, args.output)
    return True


def _parse_pipeline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="OptiBench-AGI: automated optimisation paper benchmarking",
    )
    parser.add_argument("query", nargs="?", default=None)
    parser.add_argument("--arxiv-id", default=None)
    parser.add_argument("--dataset-root", default="./dataset")
    parser.add_argument("--max-retries", type=int, default=3)
    parser.add_argument(
        "--domain",
        default=os.getenv("OPTIBENCH_DOMAIN", "continuous"),
        choices=["continuous", "all"],
    )
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument(
        "--skip-lean",
        action="store_true",
        help="Skip Lean 4 code generation/validation; still generate LaTeX prove_cot.",
    )
    parser.add_argument("--allow-empty-outline", action="store_true")
    parser.add_argument("--verbose", "-v", action="store_true")
    return parser.parse_args()


def main() -> None:
    if _dispatch_tool_subcommand():
        return

    args = interactive_pipeline_args() if len(sys.argv) == 1 else _parse_pipeline_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(name)-12s  %(levelname)-8s  %(message)s",
    )

    os.environ["OPTIBENCH_DOMAIN"] = args.domain
    from workflow import OptiBenchWorkflow

    wf = OptiBenchWorkflow(
        model=_build_model(),
        dataset_root=args.dataset_root,
        max_retries=args.max_retries,
        skip_lean=args.skip_lean,
        require_outline=not args.allow_empty_outline,
    )
    run_batch_pipeline(args, wf)


if __name__ == "__main__":
    main()
