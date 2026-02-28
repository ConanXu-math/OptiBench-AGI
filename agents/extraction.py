"""Extraction Agent – converts paper Markdown into a structured MathOutline."""

from __future__ import annotations

from typing import TYPE_CHECKING

from schema.models import MathOutline

if TYPE_CHECKING:
    from agno.models.base import Model


def create_extraction_agent(model: "Model") -> "Agent":  # noqa: F821
    """Factory that builds the Extraction Agent.

    Uses ``output_schema=MathOutline`` so the LLM response is automatically
    validated and returned as a Pydantic object.
    """
    from agno.agent import Agent

    return Agent(
        name="Extraction Agent",
        model=model,
        output_schema=MathOutline,
        instructions=[
            "You are a mathematical modelling expert.",
            "Given the Markdown text of an optimisation paper, extract the PRIMARY "
            "optimisation problem and return a MathOutline.",
            "",
            "Rules you MUST follow:",
            "1. The `objective` field must be a single, self-contained LaTeX expression "
            "(e.g. \\min_{x} f(x)).",
            "2. Each item in `constraints` must be a standalone LaTeX inequality or "
            "equality – one constraint per entry.",
            "3. `variables` should list every decision variable name (e.g. ['x', 'y']).",
            "4. `notation_table` must map EVERY symbol that appears in the objective or "
            "constraints to its dimension and meaning.",
            "5. Strictly distinguish decision variables from constants / parameters.",
            "6. All LaTeX MUST be syntactically valid and renderable – NO natural-"
            "language inside math delimiters.",
            "7. If the paper contains multiple formulations, extract the MAIN one "
            "(usually the first or most prominent).",
        ],
    )
