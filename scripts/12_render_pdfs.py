#!/usr/bin/env python3
"""
Step 12: render every report .md to PDF.
  For each results/reports/<label>/*.md (all report directories, or the labels given), converts the Markdown to
  print-styled HTML (landscape when a table has 7+ columns) and prints it with headless Chrome, Edge or Chromium
  to a .pdf beside the .md. Set $CHROME to choose the browser; under WSL the Windows Chrome / Edge builds work.
  Intermediate HTML goes to data/tmp/ and is removed afterwards. Skips cleanly when no browser is found.
"""

from __future__ import annotations

import argparse
import datetime
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.pdf_render import find_browser, render_pdf, report_markdown_files

REPO_ROOT = SCRIPT_DIR.parent
REPORTS_DIR = REPO_ROOT / "results" / "reports"


def main() -> None:
    parser = argparse.ArgumentParser(description="Render report Markdown files to PDF")
    parser.add_argument("labels", nargs="*", help="Report directories under results/reports/ (default: all)")
    args = parser.parse_args()

    browser = find_browser()
    if not browser:
        print("Step 12 skipped: no Chrome, Edge or Chromium found (set $CHROME to a browser executable).")
        return

    if args.labels:
        dirs = [REPORTS_DIR / label for label in args.labels]
        missing = [d.name for d in dirs if not d.is_dir()]
        if missing:
            print(f"ERROR: no report directory for {', '.join(missing)} under {REPORTS_DIR}", file=sys.stderr)
            sys.exit(1)
    else:
        dirs = sorted(d for d in REPORTS_DIR.iterdir() if d.is_dir()) if REPORTS_DIR.is_dir() else []

    work_dir = REPO_ROOT / "data" / "tmp"
    stamp = datetime.date.today().isoformat()
    failures = 0
    for d in dirs:
        for md in report_markdown_files(d):
            pdf = md.with_suffix(".pdf")
            footer = f"Rendered {stamp} from {d.name}/{md.name}"
            try:
                render_pdf(md, pdf, browser, work_dir, footer=footer)
                print(f"[{d.name}] {pdf.name} ({pdf.stat().st_size // 1024} KB)")
            except Exception as exc:  # one bad report must not stop the others
                failures += 1
                print(f"[{d.name}] FAILED {md.name}: {exc}", file=sys.stderr)
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
