#!/usr/bin/env python3
"""测试 LLM API 连接，诊断连不上的原因。

用法:
  python scripts/test_api_connection.py
  # 或从项目根目录:
  python -m scripts.test_api_connection
"""

from __future__ import annotations

import os
import sys
import time

# 加载 .env（与 main.py 一致）
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

# 项目根目录
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def main() -> None:
    provider = (os.getenv("OPTIBENCH_PROVIDER") or "openai").strip().lower()
    model_id = (os.getenv("OPTIBENCH_MODEL") or "gpt-4o").strip()
    base_url = (os.getenv("OPTIBENCH_BASE_URL") or "").strip()
    api_key = (os.getenv("OPTIBENCH_API_KEY") or "").strip()
    timeout_s = (os.getenv("OPTIBENCH_API_TIMEOUT") or "60").strip()

    print("=== 环境变量 ===\n")
    print(f"  OPTIBENCH_PROVIDER = {provider!r}")
    print(f"  OPTIBENCH_MODEL   = {model_id!r}")
    print(f"  OPTIBENCH_BASE_URL = {base_url!r}")
    print(f"  OPTIBENCH_API_KEY  = {(api_key[:8] + '...' + api_key[-4:]) if len(api_key) > 12 else '(未设置或过短)'}")
    print(f"  OPTIBENCH_API_TIMEOUT = {timeout_s!r}\n")

    if provider in ("openai-like", "openailike"):
        if not base_url:
            print("错误: openai-like 模式下 OPTIBENCH_BASE_URL 必须设置。")
            sys.exit(1)
        if not api_key:
            print("错误: OPTIBENCH_API_KEY 未设置。")
            sys.exit(1)
        # 检查 key 是否误带前导空格（.env 里写 OPTIBENCH_API_KEY= xxx 会带空格）
        if api_key != api_key.strip():
            print("警告: API Key 含首尾空格，已自动 strip。建议在 .env 里写成 OPTIBENCH_API_KEY=sk-xxx（等号后不要空格）。\n")
            api_key = api_key.strip()
    else:
        print("当前仅对 openai-like 做完整测试；其他 provider 可自行加逻辑。\n")

    # 1) 基础 HTTP 连通性（带重试；部分服务如 DeepSeek 对未认证请求返回 401，表示已连上）
    print("=== 1) 基础连通性（最多重试 3 次） ===\n")
    if base_url:
        import urllib.request
        from urllib.parse import urlparse
        parsed = urlparse(base_url)
        origin = f"{parsed.scheme}://{parsed.netloc}"
        last_err = None
        got_401 = False
        for attempt in range(1, 4):
            try:
                req = urllib.request.Request(origin, method="HEAD")
                urllib.request.urlopen(req, timeout=15)
                print(f"  第 {attempt} 次: 可访问 {origin}\n")
                last_err = None
                break
            except Exception as e:
                last_err = e
                err_str = str(e).lower()
                if "401" in str(e) or "unauthorized" in err_str:
                    got_401 = True
                print(f"  第 {attempt} 次: {type(e).__name__}: {e}")
                if attempt < 3:
                    time.sleep(2)
                    print("  2s 后重试…")
        if last_err is not None:
            if got_401:
                print("\n  收到 401 表示已连上服务器，但未认证或 Key 无效；请检查 OPTIBENCH_API_KEY。继续测试 Chat API。\n")
            else:
                print("\n  可能原因: 网络波动、服务端短暂不可用、或当前网络/代理与之前不同。")
                print("  若之前能用: 多跑几次本脚本或过会儿再试；换回当时用的网络/VPN 也常有帮助。\n")
    else:
        print("  未配置 BASE_URL，跳过。\n")

    # 2) 调用 Chat API
    print("=== 2) 调用 Chat API ===\n")
    if provider not in ("openai-like", "openailike"):
        print("  跳过（非 openai-like）。\n")
        return

    # 2a) 优先用 agno（与 main 一致）
    try:
        from agno.models.openai.like import OpenAILike
        use_agno = True
    except ImportError:
        use_agno = False

    if use_agno:
        kwargs = {"id": model_id, "base_url": base_url, "api_key": api_key}
        try:
            kwargs["timeout"] = float(timeout_s)
        except ValueError:
            kwargs["timeout"] = 60.0
        try:
            model = OpenAILike(**kwargs)
            response = model.response(messages=[{"role": "user", "content": "reply with one word: ok"}])
            # 兼容 agno 有时返回 dict（如 API 报错时）
            if isinstance(response, dict):
                text = (response.get("choices") or [{}])[0].get("message", {}).get("content", "") or response.get("error", response)
                text = str(text)[:200]
            else:
                text = getattr(response, "content", None) or str(response)
                text = (text or "")[:200]
            print(f"  成功. 模型回复: {text!r}\n")
        except AttributeError as e:
            if "log" in str(e).lower() or "dict" in str(e).lower():
                print("  (使用直接请求验证 Chat API)\n")
                _chat_via_urllib(base_url, model_id, api_key, timeout_s)
            else:
                _print_api_error(e)
                sys.exit(1)
        except Exception as e:
            _print_api_error(e)
            sys.exit(1)
        return

    # 2b) 无 agno 时用 urllib 发裸 POST（同样能暴露连接/SSL/401 等问题）
    print("  (未安装 agno，用 urllib 直接请求 chat/completions)\n")
    _chat_via_urllib(base_url, model_id, api_key, timeout_s)


def _chat_via_urllib(base_url: str, model_id: str, api_key: str, timeout_s: str) -> None:
    import json
    import urllib.request

    chat_url = base_url.rstrip("/") + "/chat/completions"
    body = json.dumps({
        "model": model_id,
        "messages": [{"role": "user", "content": "reply ok"}],
        "max_tokens": 10,
    }).encode("utf-8")
    req = urllib.request.Request(
        chat_url,
        data=body,
        headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {api_key}",
        },
        method="POST",
    )
    try:
        timeout = float(timeout_s) if timeout_s else 60.0
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode())
            text = (data.get("choices") or [{}])[0].get("message", {}).get("content", "")
            print(f"  成功. 模型回复: {text[:200]!r}\n")
    except Exception as e:
        _print_api_error(e)
        sys.exit(1)


def _print_api_error(e: Exception) -> None:
    print(f"  失败: {type(e).__name__}: {e}\n")
    err_lower = str(e).lower()
    if "connection" in err_lower or "timeout" in err_lower or "connect" in err_lower or "handshake" in err_lower:
        print("  可能原因: 连接/SSL 握手超时。若之前能用，多为临时波动或服务端繁忙，可稍后重试或换回当时的网络/VPN。")
    elif "401" in str(e) or "unauthorized" in err_lower or "api_key" in err_lower:
        print("  可能原因: API Key 错误或已失效；检查 .env 里 OPTIBENCH_API_KEY（等号后不要空格）。")
    elif "404" in str(e):
        print("  可能原因: base_url 或模型 path 不对，确认 OPTIBENCH_BASE_URL 是否包含完整 /v1 路径。")
    elif "429" in str(e) or "rate" in err_lower:
        print("  可能原因: 限流，稍后重试。")
    else:
        print("  建议: 根据报错检查服务端、base_url、模型名和 API Key。")
    print()


if __name__ == "__main__":
    main()
