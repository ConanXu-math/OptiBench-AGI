"""Research Agent – searches arXiv and converts PDFs to Markdown."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import arxiv

from toolkits.pdf_converter import convert_pdf_to_markdown

if TYPE_CHECKING:
    from agno.models.base import Model


def download_and_convert_paper(arxiv_id: str) -> str:
    """Download an arXiv paper by its ID and convert the PDF to Markdown.

    Args:
        arxiv_id: The arXiv identifier, e.g. ``"2301.12345"`` or
            ``"2301.12345v2"``.

    Returns:
        Markdown text of the paper, or an error message.
    """
    try:
        client = arxiv.Client()
        search = arxiv.Search(id_list=[arxiv_id])
        paper = next(client.results(search))
    except Exception as exc:
        return f"ERROR: could not fetch paper {arxiv_id} – {exc}"

    with tempfile.TemporaryDirectory(prefix="optibench_dl_") as tmpdir:
        pdf_path = paper.download_pdf(dirpath=tmpdir)
        markdown = convert_pdf_to_markdown(str(pdf_path))

    title = paper.title
    header = (
        f"# {title}\n\n"
        f"**arXiv ID:** {arxiv_id}\n\n"
        f"**Authors:** {', '.join(a.name for a in paper.authors)}\n\n"
        f"**Abstract:** {paper.summary}\n\n---\n\n"
    )
    return header + markdown


def create_research_agent(
    model: Model,
    *,
    download_dir: Path | None = None,
) -> "Agent":  # noqa: F821
    """Factory that builds the Research Agent.

    The agent has two tool sets:
    * ``ArxivTools`` for keyword search.
    * ``download_and_convert_paper`` for full-text extraction.
    """
    from agno.agent import Agent
    from agno.tools.arxiv import ArxivTools

    toolkit_kwargs = {}
    if download_dir is not None:
        toolkit_kwargs["download_dir"] = download_dir

    return Agent(
        name="Research Agent",
        model=model,
        tools=[
            ArxivTools(**toolkit_kwargs),
            download_and_convert_paper,
        ],
        instructions=[
            "You are a research assistant specialising in mathematical optimisation.",
            "When given a search query:",
            "  1. Use search_arxiv to find the most relevant papers.",
            "  2. Select the single best-matching paper.",
            "  3. Call download_and_convert_paper with the paper's arXiv ID to obtain "
            "full Markdown text.",
            "  4. Return ONLY the Markdown content of the paper (title, abstract, and "
            "the body converted from PDF). Do NOT add your own commentary.",
        ],
        markdown=True,
    )
