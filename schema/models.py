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


class SectionSelection(BaseModel):
    """Selected high-value sections for extraction."""

    selected_markdown: str = Field(
        ...,
        description="Concise markdown excerpt focusing on problem formulation/method sections",
    )
    rationale: str = Field(
        default="",
        description="Why these sections were selected",
    )


class ExtractionCritique(BaseModel):
    """Critic output for extracted MathOutline."""

    score: int = Field(
        ...,
        ge=0,
        le=100,
        description="Quality score of extracted outline",
    )
    is_acceptable: bool = Field(
        ...,
        description="Whether extraction quality is acceptable for downstream coding",
    )
    issues: list[str] = Field(
        default_factory=list,
        description="Specific issues found in objective/constraints/variables",
    )
    fix_prompt: str = Field(
        default="",
        description="Actionable prompt to guide next extraction attempt",
    )


class QualityGateResult(BaseModel):
    """Overall quality-gate result before final save."""

    passed: bool = Field(..., description="Whether benchmark item passes quality gate")
    checks: dict[str, str] = Field(
        default_factory=dict,
        description="Per-check status map: ok/warn/fail",
    )
    notes: list[str] = Field(
        default_factory=list,
        description="Human-readable notes for failures/warnings",
    )
