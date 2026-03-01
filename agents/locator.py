"""Locator Agent – finds high-value sections for math extraction."""

from __future__ import annotations

from typing import TYPE_CHECKING

from schema.models import SectionSelection

if TYPE_CHECKING:
    from agno.models.base import Model


def create_locator_agent(model: "Model") -> "Agent":  # noqa: F821
    """Factory that builds the section locator agent."""
    from agno.agent import Agent

    return Agent(
        name="Locator Agent",
        model=model,
        output_schema=SectionSelection,
        instructions=[
            "You are a paper structure analyst for optimization papers.",
            "Given full markdown, select ONLY high-value sections for math extraction:",
            "- abstract",
            "- problem formulation",
            "- objective/constraint definitions",
            "- algorithm statements with equations",
            "- convergence theorem statements (without long proofs)",
            "Return concise markdown excerpt in `selected_markdown`.",
            "Do not include references, appendices, or unrelated background.",
            "Target 4k-12k characters unless source is shorter.",
        ],
    )

