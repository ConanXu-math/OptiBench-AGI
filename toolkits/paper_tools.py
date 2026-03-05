"""Paper 相关工具：arXiv 搜索、下载、PDF 转 Markdown。

可在 CLI 子命令或 Agent 中复用。
支持环境变量 OPTIBENCH_ARXIV_TIMEOUT、OPTIBENCH_ARXIV_BASE_URL（镜像）。
"""

from __future__ import annotations

import os
from pathlib import Path
from urllib.parse import urlparse
from typing import Any

import arxiv

from toolkits.pdf_converter import convert_pdf_to_markdown

# 默认 API 端点（可被 OPTIBENCH_ARXIV_BASE_URL 覆盖为镜像）
_DEFAULT_ARXIV_API_BASE = "https://export.arxiv.org"


def _arxiv_client() -> arxiv.Client:
    """创建已配置超时与可选镜像的 arXiv Client。"""
    base = (os.getenv("OPTIBENCH_ARXIV_BASE_URL") or "").strip().rstrip("/")
    if base:
        arxiv.Client.query_url_format = f"{base}/api/query?{{}}"
    else:
        arxiv.Client.query_url_format = f"{_DEFAULT_ARXIV_API_BASE}/api/query?{{}}"

    timeout_s = (os.getenv("OPTIBENCH_ARXIV_TIMEOUT") or "120").strip()
    try:
        timeout_f = float(timeout_s)
    except ValueError:
        timeout_f = 120.0

    client = arxiv.Client()
    orig_get = client._session.get

    def _get_with_timeout(*args: Any, **kwargs: Any) -> Any:
        kwargs.setdefault("timeout", timeout_f)
        return orig_get(*args, **kwargs)

    client._session.get = _get_with_timeout
    return client


def _arxiv_download_domain() -> str:
    """用于 PDF 下载的域名（Result.download_pdf(download_domain=...)）。"""
    base = (os.getenv("OPTIBENCH_ARXIV_BASE_URL") or _DEFAULT_ARXIV_API_BASE).strip()
    if not base.startswith(("http://", "https://")):
        base = "https://" + base
    parsed = urlparse(base)
    return parsed.netloc or "export.arxiv.org"


def get_arxiv_info(arxiv_id: str) -> dict[str, Any] | None:
    """根据 arXiv ID 获取论文元数据（标题、作者、摘要等），不下载。

    Returns:
        包含 title, arxiv_id, summary, pdf_url, authors 的字典；未找到则返回 None。
    """
    client = _arxiv_client()
    search = arxiv.Search(id_list=[arxiv_id])
    try:
        paper = next(client.results(search))
    except StopIteration:
        return None
    raw_id = getattr(paper, "entry_id", None) or ""
    aid = raw_id.split("/abs/")[-1].strip() if isinstance(raw_id, str) else arxiv_id
    return {
        "title": paper.title,
        "arxiv_id": aid,
        "summary": paper.summary or "",
        "pdf_url": getattr(paper, "pdf_url", None) or f"https://arxiv.org/pdf/{aid}.pdf",
        "authors": [a.name for a in paper.authors],
        "published": getattr(paper, "published", None),
    }


def search_arxiv(
    query: str,
    max_results: int = 10,
    *,
    domain: str | None = None,
) -> list[dict[str, Any]]:
    """在 arXiv 上按关键词搜索论文。

    Args:
        query: 搜索关键词，如 "convex optimization"。
        max_results: 最多返回条数。

    Returns:
        列表，每项为 {"title", "arxiv_id", "summary", "pdf_url", "authors"}。
    """
    if domain is None:
        domain = os.getenv("OPTIBENCH_DOMAIN", "continuous")
    domain = (domain or "all").strip().lower()

    # Pull a larger pool first, then filter by domain to keep enough candidates.
    candidate_pool = max(max_results * 5, max_results)
    client = _arxiv_client()
    search = arxiv.Search(query=query, max_results=candidate_pool)
    out: list[dict[str, Any]] = []
    for p in client.results(search):
        categories = list(getattr(p, "categories", []) or [])
        if domain in {"continuous", "continuous-optimization"} and not _is_continuous_optimization_paper(
            title=p.title,
            summary=p.summary or "",
            categories=categories,
        ):
            continue
        # entry_id 形如 https://arxiv.org/abs/1406.0899v4
        raw_id = getattr(p, "entry_id", None) or getattr(p, "arxiv_id", str(p))
        arxiv_id = raw_id.split("/abs/")[-1].strip() if isinstance(raw_id, str) else str(raw_id)
        out.append({
            "title": p.title,
            "arxiv_id": arxiv_id,
            "summary": (p.summary or "")[:500],
            "pdf_url": p.pdf_url,
            "authors": [a.name for a in p.authors],
            "categories": categories,
        })
        if len(out) >= max_results:
            break
    return out


def _is_continuous_optimization_paper(
    *,
    title: str,
    summary: str,
    categories: list[str],
) -> bool:
    """Heuristic filter for traditional continuous optimization papers."""
    text = f"{title}\n{summary}".lower()
    cats = {c.lower() for c in categories}

    # Prefer classical optimization categories.
    category_match = bool(cats.intersection({"math.oc", "math.na", "math.ap"}))

    include_terms = [
        "convex optimization",
        "convex optimisation",
        "non-convex optimization",
        "nonconvex optimization",
        "lagrangian",
        "kkt",
        "duality",
        "interior-point",
        "trust-region",
        "line search",
        "gradient method",
        "first-order method",
        "second-order",
        "variational inequality",
        "semidefinite",
        "quadratic programming",
        "linear programming",
    ]
    include_match = any(t in text for t in include_terms)

    # Exclude mostly ML-centric papers that often dominate naive search.
    exclude_terms = [
        "reinforcement learning",
        "deep learning",
        "neural network",
        "transformer",
        "diffusion model",
        "language model",
        "llm",
        "bert",
        "gpt",
    ]
    exclude_match = any(t in text for t in exclude_terms)

    return (category_match or include_match) and not exclude_match


def download_arxiv_pdf(arxiv_id: str, out_dir: str | Path) -> Path:
    """按 arXiv ID 下载 PDF 到指定目录。

    Args:
        arxiv_id: 如 "1406.0899" 或 "1406.0899v4"。
        out_dir: 保存目录，若不存在会创建。

    Returns:
        下载后的 PDF 文件路径。

    Raises:
        ValueError: 未找到该 ID 的论文。
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    client = _arxiv_client()
    search = arxiv.Search(id_list=[arxiv_id])
    try:
        paper = next(client.results(search))
    except StopIteration:
        raise ValueError(f"No paper found for arXiv ID: {arxiv_id}") from None
    pdf_path = paper.download_pdf(dirpath=str(out_dir), download_domain=_arxiv_download_domain())
    return Path(pdf_path)


def download_paper(
    arxiv_id: str,
    out_dir: str | Path,
    *,
    convert_to_md: bool = True,
    save_md_file: bool = True,
) -> dict[str, Any]:
    """下载论文并可选择转为 Markdown。

    Args:
        arxiv_id: arXiv 标识。
        out_dir: 输出目录，会创建 arxiv_id 子目录。
        convert_to_md: 是否用 marker/nougat 转为 Markdown。
        save_md_file: 转为 Markdown 时是否写入 .md 文件（否则只返回内容）。

    Returns:
        {"pdf_path", "md_path" | None, "md_content" | None, "title", "authors"}。
    """
    out_dir = Path(out_dir)
    paper_dir = out_dir / arxiv_id.replace("/", "_")
    paper_dir.mkdir(parents=True, exist_ok=True)

    client = _arxiv_client()
    search = arxiv.Search(id_list=[arxiv_id])
    try:
        paper = next(client.results(search))
    except StopIteration:
        raise ValueError(f"No paper found for arXiv ID: {arxiv_id}") from None

    pdf_path = paper.download_pdf(dirpath=str(paper_dir), download_domain=_arxiv_download_domain())
    pdf_path = Path(pdf_path)
    result: dict[str, Any] = {
        "pdf_path": pdf_path,
        "md_path": None,
        "md_content": None,
        "title": paper.title,
        "authors": [a.name for a in paper.authors],
        "summary": paper.summary,
    }

    if convert_to_md:
        md_content = convert_pdf_to_markdown(str(pdf_path))
        header = (
            f"# {paper.title}\n\n"
            f"**arXiv ID:** {arxiv_id}\n\n"
            f"**Authors:** {', '.join(a.name for a in paper.authors)}\n\n"
            f"**Abstract:** {paper.summary}\n\n---\n\n"
        )
        result["md_content"] = header + md_content
        if save_md_file:
            md_path = paper_dir / "paper.md"
            md_path.write_text(result["md_content"], encoding="utf-8")
            result["md_path"] = md_path

    return result
