"""Coding Agent – generates Python + pseudocode from a MathOutline."""

from __future__ import annotations

from typing import TYPE_CHECKING

from schema.models import CodingOutput

if TYPE_CHECKING:
    from agno.models.base import Model


def create_coding_agent(model: "Model") -> "Agent":  # noqa: F821
    """Factory that builds the Coding Agent.

    Returns structured ``CodingOutput`` (pseudocode + pycode).
    """
    from agno.agent import Agent

    return Agent(
        name="Coding Agent",
        model=model,
        output_schema=CodingOutput,
        instructions=[
            "You are an algorithm engineer who translates mathematical optimisation "
            "formulations into executable Python code.",
            "",
            "Given a MathOutline (JSON), produce TWO artefacts:",
            "",
            "1. **pseudocode** – concise, language-agnostic pseudocode of the algorithm.",
            "2. **pycode** – a COMPLETE, self-contained Python script that satisfies ALL "
            "of the following hard constraints:",
            "",
            "Hard constraints for `pycode`:",
            "  a) Import only `numpy` and `scipy.optimize` (plus stdlib).",
            "  b) Define a function `generate_synthetic_data(seed: int = 42) -> dict` "
            "that creates a feasible test instance honouring every constraint in the "
            "MathOutline.  The returned dict maps variable/parameter names to numpy "
            "arrays.",
            "  c) Define a function `solve(data: dict) -> dict` that solves the "
            "optimisation problem using scipy.optimize (e.g. minimize, linprog, "
            "milp) and returns a dict with at least 'optimal_value' and 'solution'.",
            "  d) Include an `if __name__ == '__main__':` block that calls "
            "generate_synthetic_data(), then solve(), and prints the results.",
            "  e) The script must run without errors on a clean Python 3.13 environment "
            "with numpy and scipy installed.",
            "  f) Do NOT use any plotting or GUI libraries.",
            "  g) Do NOT read external files; all data comes from generate_synthetic_data.",
        ],
    )
