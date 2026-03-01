"""Quality-gate checks before saving benchmark artifacts."""

from __future__ import annotations

from schema.models import BenchmarkItem, QualityGateResult


def evaluate_quality_gate(item: BenchmarkItem, *, skip_lean: bool) -> QualityGateResult:
    checks: dict[str, str] = {}
    notes: list[str] = []

    # Extraction quality
    if item.outline.objective.strip() and (
        item.outline.variables or item.outline.constraints
    ):
        checks["outline"] = "ok"
    else:
        checks["outline"] = "fail"
        notes.append("Outline missing objective or variables/constraints.")

    # Coding quality
    if item.pycode.strip() and item.pseudocode.strip():
        checks["coding"] = "ok"
    else:
        checks["coding"] = "fail"
        notes.append("Coding artifacts missing pycode/pseudocode.")

    # Prove quality
    if item.prove_cot.strip():
        checks["prove_cot"] = "ok"
    else:
        checks["prove_cot"] = "fail"
        notes.append("prove_cot is empty.")

    # Lean quality (optional)
    if skip_lean:
        checks["lean4"] = "skipped"
    else:
        if item.lean4_formal.strip():
            checks["lean4"] = "ok"
        else:
            checks["lean4"] = "fail"
            notes.append("lean4_formal is empty while skip_lean=False.")

    passed = not any(v == "fail" for v in checks.values())
    return QualityGateResult(passed=passed, checks=checks, notes=notes)

