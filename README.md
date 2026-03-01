# OptiBench-AGI

Automated optimisation-paper benchmarking pipeline built on the **Agno** multi-agent framework.

The system searches arXiv, extracts mathematical formulations, generates validated Python solvers, and produces Lean 4 formal statements — all orchestrated by four specialised AI agents.

## Architecture

```
┌──────────────┐     ┌──────────────────┐     ┌──────────────┐     ┌────────────────────┐
│ Research Agent│ ──▶ │ Extraction Agent  │ ──▶ │ Coding Agent │ ──▶ │ Formalization Agent │
│ (arXiv+PDF)  │     │ (→ MathOutline)   │     │ (Python+pseudo)   │ (Lean 4)            │
└──────────────┘     └──────────────────┘     └───────┬──────┘     └─────────┬──────────┘
                                                      │ retry ≤ 3           │ retry ≤ 3
                                                      ▼                     ▼
                                               validate_python        validate_lean
```

## Quick Start

```bash
# 1. Install (core dependencies)
uv pip install -e .

# 2. Install PDF converter (choose one)
uv pip install -e ".[marker]"   # marker-pdf (recommended)
# or
uv pip install -e ".[nougat]"   # nougat-ocr

# 3. Set your API key
export OPENAI_API_KEY="sk-..."

# 4. Run
python main.py "convex optimisation relaxation techniques"

# Or target a specific paper
python main.py --arxiv-id 2301.12345
```

## Tools (no API key required)

Standalone commands for searching and downloading papers, or converting PDFs:

```bash
# Search arXiv and print results
python main.py search "convex optimisation" --max 5

# Download a paper by arXiv ID (PDF + optional Markdown)
python main.py download 1406.0899
python main.py download 1406.0899 --output ./papers --no-convert   # PDF only

# Convert a local PDF to Markdown (requires marker or nougat)
python main.py convert-pdf path/to/paper.pdf -o paper.md

# List processed papers in dataset/
python main.py list
python main.py list --dataset-root ./dataset

# Show metadata for an arXiv ID (title, authors, abstract)
python main.py info 1406.0899
```

## Configuration

All configuration is via environment variables:

| Variable | Default | Description |
|---|---|---|
| `OPENAI_API_KEY` | — | OpenAI API key (required for default provider) |
| `OPTIBENCH_PROVIDER` | `openai` | LLM provider: `openai`, `anthropic`, `google` |
| `OPTIBENCH_MODEL` | `gpt-4o` | Model identifier |
| `OPTIBENCH_PY_TIMEOUT` | `120` | Python validation timeout (seconds) |
| `OPTIBENCH_LEAN_TIMEOUT` | `300` | Lean 4 build timeout (seconds) |
| `OPTIBENCH_LEAN_PROJECT_DIR` | — | Pre-built Lean+Mathlib project path (speeds up validation) |
| `OPTIBENCH_MARKER_TIMEOUT` | `300` | marker PDF conversion timeout |

## Output Structure

Each processed paper produces a folder under `./dataset/<arxiv_id>/`:

```
dataset/
└── 2301.12345/
    ├── benchmark.json   # Complete BenchmarkItem
    ├── outline.json     # MathOutline only
    ├── solve.py         # Validated Python solver
    ├── pseudocode.txt   # Algorithm pseudocode
    └── Formal.lean      # Lean 4 theorem statements
```

## Data Models

- **NotationItem** — `(symbol, dimension, description)`
- **MathOutline** — `(objective, constraints, variables, notation_table)`
- **BenchmarkItem** — `(paper_name, arxiv_id, outline, prove_cot, pseudocode, pycode, lean4_formal)`

## Lean 4 Validation

For Lean 4 validation you need:

1. **elan** — the Lean toolchain installer
2. A project with **Mathlib** dependency

Fastest approach: pre-build a Lean project once and set `OPTIBENCH_LEAN_PROJECT_DIR`:

```bash
mkdir lean-check && cd lean-check
lake init OptiBenchCheck math
lake build   # first build downloads Mathlib (~10 min)
export OPTIBENCH_LEAN_PROJECT_DIR=$(pwd)
```

## Project Layout

```
├── main.py              # CLI entry point
├── workflow.py           # OptiBenchWorkflow orchestrator
├── schema/
│   └── models.py         # Pydantic data models
├── agents/
│   ├── research.py       # arXiv search + PDF conversion
│   ├── extraction.py     # Markdown → MathOutline
│   ├── coding.py         # MathOutline → Python + pseudocode
│   └── formalization.py  # MathOutline → Lean 4
├── toolkits/
│   ├── validators.py     # validate_python_code / validate_lean_code
│   └── pdf_converter.py  # marker / nougat PDF→Markdown
└── pyproject.toml
```
