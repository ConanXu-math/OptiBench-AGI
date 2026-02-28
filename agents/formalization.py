"""Formalization Agent – generates Lean 4 theorem statements."""

from __future__ import annotations

from typing import TYPE_CHECKING

from schema.models import FormalizationOutput

if TYPE_CHECKING:
    from agno.models.base import Model


def create_formalization_agent(model: "Model") -> "Agent":  # noqa: F821
    """Factory that builds the Formalization Agent.

    Returns structured ``FormalizationOutput`` (prove_cot + lean4_formal).
    """
    from agno.agent import Agent

    return Agent(
        name="Formalization Agent",
        model=model,
        output_schema=FormalizationOutput,
        instructions=[
            "You are a Lean 4 formalisation expert with deep knowledge of Mathlib.",
            "",
            "Given a MathOutline (JSON), produce:",
            "1. `prove_cot` – a chain-of-thought explaining your formalisation strategy.",
            "2. `lean4_formal` – valid Lean 4 source code satisfying ALL rules below.",
            "",
            "Rules for `lean4_formal`:",
            "  a) Start with: import Mathlib.Analysis.Convex.Function",
            "  b) Open relevant namespaces as needed (e.g. `open scoped NNReal`).",
            "  c) Define the ambient space as `EuclideanSpace ℝ (Fin n)` where n is a "
            "universe variable or explicit Nat.",
            "  d) State each theorem / definition precisely. Proofs MUST use `sorry` – "
            "we only care about well-typed *statements*.",
            "  e) The code must compile with `lake build` in a project that depends on "
            "Mathlib (commit compatible with lean4 v4.14.0).",
            "  f) Do NOT invent Mathlib names that do not exist; stick to the real API.",
            "  g) Keep the file self-contained (single .lean file).",
        ],
    )
