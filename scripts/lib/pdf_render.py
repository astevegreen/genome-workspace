"""Markdown report -> print-styled HTML -> PDF (headless Chrome / Edge / Chromium), used by Step 12."""

from __future__ import annotations

import html
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import markdown

# A report whose widest table has at least this many columns prints in landscape (the hand-made 22 Nov PDFs
# used landscape for the 9-column Step 10 and couple reports and portrait for the 6-column card and summary).
LANDSCAPE_MIN_COLUMNS = 7

_TABLE_ROW = re.compile(r"^\s*\|")
_TABLE_SEP = re.compile(r"^\s*\|?\s*:?-{3,}:?\s*(\|\s*:?-{3,}:?\s*)*\|?\s*$")
_LIST_ITEM = re.compile(r"^\s*([-*+]|\d+[.)])\s+")
_FENCE = re.compile(r"^\s*(```|~~~)")

# Browsers tried in order when $CHROME is unset: Linux builds first, then Windows builds through WSL interop.
BROWSER_CANDIDATES = (
    "chromium",
    "chromium-browser",
    "google-chrome",
    "google-chrome-stable",
    "microsoft-edge",
    "/mnt/c/Program Files/Google/Chrome/Application/chrome.exe",
    "/mnt/c/Program Files (x86)/Google/Chrome/Application/chrome.exe",
    "/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe",
    "/mnt/c/Program Files/Microsoft/Edge/Application/msedge.exe",
)


def normalize_markdown(text: str) -> str:
    """Insert the blank line Python-Markdown needs before a list or table that follows a paragraph line.

    The pipeline writes GitHub-style Markdown, where "**Label:**" directly followed by "- item" or "| a | b |"
    still starts a list or table; Python-Markdown would fold those lines into the paragraph instead.
    Fenced code blocks are left untouched.
    """
    out: list[str] = []
    in_fence = False
    prev = ""
    for line in text.splitlines():
        if _FENCE.match(line):
            if not in_fence and prev.strip():
                out.append("")
            in_fence = not in_fence
            out.append(line)
            prev = line
            continue
        if not in_fence and prev.strip() and not prev.startswith((" ", "\t")):
            starts_list = bool(_LIST_ITEM.match(line)) and not _LIST_ITEM.match(prev)
            starts_table = bool(_TABLE_ROW.match(line)) and not _TABLE_ROW.match(prev)
            if (starts_list or starts_table) and not prev.lstrip().startswith(">"):
                out.append("")
        out.append(line if in_fence else escape_star_alleles(line))
        prev = line
    return "\n".join(out) + "\n"


_STAR_ALLELE = re.compile(r"(?<![*\\])\*(?=\d)")


def escape_star_alleles(line: str) -> str:
    """Escape the '*' of star-allele names (X*99:01, *9/*99) outside code spans so they are not read as italics.

    A '*' directly before a digit is never emphasis in the pipeline's reports; '**' (bold) is left alone.
    """
    parts = line.split("`")
    if len(parts) % 2 == 0:  # unbalanced backticks: leave the line as written
        return line
    return "`".join(_STAR_ALLELE.sub(r"\\*", p) if i % 2 == 0 else p for i, p in enumerate(parts))


def max_table_columns(text: str) -> int:
    """Widest Markdown table in the text (columns counted on each table's separator row)."""
    widest = 0
    for line in text.splitlines():
        if _TABLE_SEP.match(line) and "|" in line:
            cells = line.strip().strip("|").split("|")
            widest = max(widest, len(cells))
    return widest


def page_orientation(text: str) -> str:
    return "landscape" if max_table_columns(text) >= LANDSCAPE_MIN_COLUMNS else "portrait"


def markdown_to_html_body(text: str) -> str:
    return markdown.markdown(normalize_markdown(text), extensions=["tables", "fenced_code", "sane_lists"])


_CSS = """
  @page { size: letter %(orientation)s; margin: 0.55in 0.55in 0.6in 0.55in; }
  body { font-family: "Segoe UI", Roboto, Helvetica, Arial, sans-serif; font-size: %(font)s; line-height: 1.42;
         color: #1f2937; margin: 0; }
  h1 { font-size: 16pt; color: #0f172a; border-bottom: 2px solid #2563eb; padding-bottom: 5px; margin: 0 0 10px; }
  h2 { font-size: 12pt; color: #1e3a8a; border-bottom: 1px solid #cbd5e1; padding-bottom: 3px; margin: 16px 0 8px; }
  h3 { font-size: 10.2pt; color: #334155; margin: 12px 0 5px; }
  h4 { font-size: 9.6pt; color: #334155; margin: 10px 0 4px; }
  h1, h2, h3, h4 { break-after: avoid; }
  p { margin: 5px 0 8px; }
  ul, ol { margin: 4px 0 8px; padding-left: 20px; }
  li { margin-bottom: 3px; }
  blockquote { margin: 8px 0 12px; padding: 6px 12px; background: #eff6ff; border-left: 4px solid #2563eb;
               color: #1e3a8a; break-inside: avoid; }
  blockquote p { margin: 3px 0; }
  code { font-family: Consolas, "Courier New", monospace; font-size: 0.9em; background: #f1f5f9; padding: 1px 4px;
         border-radius: 3px; color: #0f172a; overflow-wrap: break-word; }
  pre { background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 4px; padding: 8px 10px; font-size: 8pt;
        white-space: pre-wrap; break-inside: avoid; }
  pre code { background: none; padding: 0; }
  hr { border: none; border-top: 1px solid #e2e8f0; margin: 12px 0; }
  table { width: 100%%; border-collapse: collapse; margin: 6px 0 12px; font-size: %(table_font)s; }
  thead { display: table-header-group; }
  tr { break-inside: avoid; }
  th, td { border: 1px solid #cbd5e1; padding: 3px 5px; vertical-align: top; overflow-wrap: break-word; }
  th { background: #f1f5f9; color: #0f172a; font-weight: 600; }
  tbody tr:nth-child(even) { background: #f8fafc; }
  .footer { margin-top: 18px; padding-top: 4px; border-top: 1px solid #e2e8f0; color: #64748b; font-size: 7.5pt; }
"""


def build_html_document(title: str, text: str, footer: str = "") -> str:
    """Full HTML page for one Markdown report, with print CSS sized for its orientation."""
    orientation = page_orientation(text)
    css = _CSS % {
        "orientation": orientation,
        "font": "8.5pt" if orientation == "landscape" else "9.5pt",
        "table_font": "7.4pt" if orientation == "landscape" else "8.2pt",
    }
    foot = f'<div class="footer">{html.escape(footer)}</div>' if footer else ""
    return (
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
        f"<title>{html.escape(title)}</title>\n<style>{css}</style>\n</head>\n<body>\n"
        f"{markdown_to_html_body(text)}\n{foot}\n</body>\n</html>\n"
    )


def report_markdown_files(report_dir: Path) -> list[Path]:
    """The report .md files directly in one results/reports/<label>/ directory (PharmCAT's own output is HTML)."""
    return sorted(p for p in report_dir.glob("*.md") if p.is_file())


def find_browser(env: dict[str, str] | None = None) -> str | None:
    """$CHROME if set, else the first Chromium-family browser found (Linux, then Windows through WSL)."""
    env = os.environ if env is None else env
    if env.get("CHROME"):
        return env["CHROME"]
    for cand in BROWSER_CANDIDATES:
        if cand.startswith("/"):
            if Path(cand).exists():
                return cand
        elif shutil.which(cand):
            return shutil.which(cand)
    return None


def _is_windows_exe(browser: str) -> bool:
    return browser.lower().endswith(".exe")


def _browser_path(path: Path, browser: str) -> str:
    """Path as the browser sees it: Windows builds run through WSL interop need a Windows path."""
    if _is_windows_exe(browser) and shutil.which("wslpath"):
        return subprocess.run(["wslpath", "-w", str(path)], capture_output=True, text=True, check=True).stdout.strip()
    return str(path)


def _file_url(path: Path, browser: str) -> str:
    p = _browser_path(path, browser)
    if _is_windows_exe(browser):
        return "file:///" + p.replace("\\", "/").lstrip("/")
    return path.resolve().as_uri()


def chrome_print_command(browser: str, html_path: Path, pdf_path: Path, profile_dir: Path) -> list[str]:
    return [
        browser,
        "--headless",
        "--disable-gpu",
        "--no-first-run",
        "--no-default-browser-check",
        "--disable-extensions",
        "--no-pdf-header-footer",
        f"--user-data-dir={_browser_path(profile_dir, browser)}",
        f"--print-to-pdf={_browser_path(pdf_path, browser)}",
        _file_url(html_path, browser),
    ]


def render_pdf(md_path: Path, pdf_path: Path, browser: str, work_dir: Path, footer: str = "", timeout: int = 120) -> None:
    """Render one Markdown report to PDF. The PDF is written to a temp file and moved into place only when valid."""
    text = md_path.read_text(encoding="utf-8")
    work_dir.mkdir(parents=True, exist_ok=True)
    # a Windows browser can hold its profile files open for a moment after exiting, so cleanup must not fail the run
    with tempfile.TemporaryDirectory(dir=work_dir, prefix="pdf_", ignore_cleanup_errors=True) as tmp:
        tmp_dir = Path(tmp)
        html_path = tmp_dir / f"{md_path.stem}.html"
        with open(html_path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(build_html_document(md_path.stem, text, footer))
        tmp_pdf = tmp_dir / f"{md_path.stem}.pdf"
        cmd = chrome_print_command(browser, html_path, tmp_pdf, tmp_dir / "profile")
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if not tmp_pdf.exists() or tmp_pdf.stat().st_size == 0 or tmp_pdf.read_bytes()[:5] != b"%PDF-":
            tail = (proc.stderr or proc.stdout or f"exit {proc.returncode}").strip()[-800:]
            raise RuntimeError(f"browser did not write a PDF for {md_path.name}: {tail}")
        shutil.move(str(tmp_pdf), str(pdf_path))
