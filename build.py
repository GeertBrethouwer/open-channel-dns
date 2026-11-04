#!/usr/bin/env python3
"""Build the Open Channel DNS web app.

Sources:  src/app.html (page + UI) and src/solver.js (the DNS solver).
Outputs:
  docs/index.html        standalone page with its own <html>/<head>; served by GitHub Pages
  open-channel-dns.html  page body only, for publishing as a Claude artifact
                         (the artifact host adds the document skeleton itself)

Usage:  python3 build.py
"""
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent
SOLVER_TAG = "/*__SOLVER__*/"
EXPORT_LINE = 'if (typeof module !== "undefined") module.exports = { OpenChannelDNS, FFT1 };'

# Head and base styles for the standalone page (what the artifact host would otherwise provide).
HEAD = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<meta name="description" content="Interactive direct numerical simulation of turbulent open channel flow, running in the browser.">
<style>
:root { padding-top: env(safe-area-inset-top, 0px); padding-bottom: env(safe-area-inset-bottom, 0px); }
body { margin: 0; }
img { max-width: 100%; }
[hidden] { display: none !important; }
</style>
"""


def page_body() -> str:
    app = (ROOT / "src/app.html").read_text(encoding="utf-8")
    solver = (ROOT / "src/solver.js").read_text(encoding="utf-8").replace(EXPORT_LINE, "")
    if SOLVER_TAG not in app:
        raise SystemExit("src/app.html has no " + SOLVER_TAG + " placeholder")
    if "</script" in solver:
        raise SystemExit("src/solver.js must not contain '</script'")
    # Drop head-only tags kept in app.html for local testing; each build adds its own.
    lines = [
        l for l in app.replace(SOLVER_TAG, solver).split("\n")
        if not l.startswith('<meta charset') and not l.startswith('<meta name="viewport"')
    ]
    return "\n".join(lines)


def main() -> None:
    body = page_body()

    # Claude artifact version: content only.
    (ROOT / "open-channel-dns.html").write_text(body, encoding="utf-8")

    # Standalone version for GitHub Pages: move <title>, <link> and <style> into <head>.
    head_part, rest = body.split("<header", 1)
    head_part = "\n".join(
        l for l in head_part.split("\n") if not l.startswith('<meta name="claude:')
    )
    html = HEAD + head_part.strip() + "\n</head>\n<body>\n<header" + rest + "\n</body>\n</html>\n"
    docs = ROOT / "docs"
    docs.mkdir(exist_ok=True)
    (docs / "index.html").write_text(html, encoding="utf-8")
    (docs / ".nojekyll").write_text("", encoding="utf-8")

    print(f"built docs/index.html ({len(html):,} bytes) and open-channel-dns.html ({len(body):,} bytes)")


if __name__ == "__main__":
    main()
