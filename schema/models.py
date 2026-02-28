"""Pydantic data models for the OptiBench-AGI pipeline."""

from __future__ import annotations

from pydantic import BaseModel, Field


class NotationItem(BaseModel):
    """A single entry in the notation table."""

    symbol: str = Field(..., description="LaTeX representation, e.g. \\mathbf{x}")
    dimension: str = Field(..., description="Type / dimension, e.g. \\mathbb{R}^n")
    description: str = Field(..., description="Plain-text meaning of the symbol")


class MathOutline(BaseModel):
    """Structured mathematical formulation extracted from a paper."""

    objective: str = Field(..., description="Objective function in LaTeX")
    constraints: list[str] = Field(
        default_factory=list,
        description="Each constraint as a standalone LaTeX expression",
    )
    variables: list[str] = Field(
        default_factory=list,
        description="Decision variable names used in the formulation",
    )
    notation_table: list[NotationItem] = Field(
        default_factory=list,
        description="Complete notation table mapping symbols to meanings",
    )


class CodingOutput(BaseModel):
    """Output schema for the Coding Agent."""

    pseudocode: str = Field(..., description="Algorithm pseudocode in plain text")
    pycode: str = Field(
        ...,
        description=(
            "Complete, self-contained Python script using numpy/scipy. "
            "Must include generate_synthetic_data() and solve() functions."
        ),
    )


class FormalizationOutput(BaseModel):
    """Output schema for the Formalization Agent."""

    prove_cot: str = Field(
        ..., description="Chain-of-thought reasoning explaining the formalization"
    )
    lean4_formal: str = Field(
        ...,
        description=(
            "Lean 4 source code that imports Mathlib, "
            "states the theorem, and uses sorry for proofs"
        ),
    )


class BenchmarkItem(BaseModel):
    """Complete benchmark record for one optimization paper."""

    paper_name: str = Field(..., description="Title of the source paper")
    arxiv_id: str = Field(default="", description="ArXiv paper identifier")
    outline: MathOutline = Field(..., description="Extracted mathematical formulation")
    prove_cot: str = Field(default="", description="Chain-of-thought proof reasoning")
    pseudocode: str = Field(default="", description="Algorithm pseudocode")
    pycode: str = Field(default="", description="Validated Python implementation")
    lean4_formal: str = Field(default="", description="Lean 4 formal statement")
