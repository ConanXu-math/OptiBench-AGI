"""Extraction Critic Agent – evaluates MathOutline quality."""

from __future__ import annotations

from typing import TYPE_CHECKING

from schema.models import ExtractionCritique

if TYPE_CHECKING:
    from agno.models.base import Model


def create_extraction_critic_agent(model: "Model") -> "Agent":  # noqa: F821
    """Factory that builds the extraction critic agent."""
    from agno.agent import Agent

    return Agent(
        name="Extraction Critic Agent",
        model=model,
        output_schema=ExtractionCritique,
        instructions=[
            "You are a strict reviewer for optimization-model extraction quality.",
            "Given source markdown excerpt and extracted MathOutline:",
            "1) score quality from 0-100",
            "2) decide `is_acceptable`",
            "3) list concrete issues",
            "4) provide a concise `fix_prompt` for next extraction iteration",
            "Reject outputs with missing objective, empty variables, invalid/ambiguous constraints,",
            "or obvious mismatch to source equations.",
        ],
    )

