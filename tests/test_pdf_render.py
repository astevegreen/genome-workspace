"""Unit tests for Step 12 Markdown -> HTML -> PDF rendering (synthetic reports only)."""

from __future__ import annotations

import stat
import sys
from pathlib import Path

import pytest

SCRIPT_DIR = Path(__file__).resolve().parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lib.pdf_render import (
    build_html_document,
    chrome_print_command,
    escape_star_alleles,
    find_browser,
    markdown_to_html_body,
    max_table_columns,
    normalize_markdown,
    page_orientation,
    render_pdf,
    report_markdown_files,
)


def _table(cols: int, rows: int = 1) -> str:
    head = "| " + " | ".join(f"H{i}" for i in range(cols)) + " |"
    sep = "|" + "|".join("---" for _ in range(cols)) + "|"
    body = ["| " + " | ".join(f"r{r}c{i}" for i in range(cols)) + " |" for r in range(rows)]
    return "\n".join([head, sep, *body])


def test_normalize_inserts_blank_line_before_list_and_table_after_paragraph():
    text = "**Findings:**\n- one\n- two\nIntro line\n" + _table(2)
    out = normalize_markdown(text).splitlines()
    assert out[:4] == ["**Findings:**", "", "- one", "- two"]
    assert out[4:6] == ["Intro line", ""]
    assert out[6].startswith("| H0")


def test_normalize_leaves_consecutive_items_fences_and_quotes_alone():
    text = "- a\n- b\n\n```\nplain\n- not a list\n| not | a table |\n```\n> quote\n> - quoted item\n"
    out = normalize_markdown(text)
    assert "- a\n- b\n" in out
    assert "plain\n- not a list\n| not | a table |\n```" in out
    assert "> quote\n> - quoted item" in out


def test_lists_and_tables_render_after_paragraph_line():
    body = markdown_to_html_body("**Label:**\n- first\n- second\n\nText\n" + _table(3, 2))
    assert "<ul>" in body and body.count("<li>") == 2
    assert "<table>" in body and body.count("<tr>") == 3


def test_escape_star_alleles_outside_code_only():
    assert escape_star_alleles("Gene X*9/*99 and Y*77:01") == r"Gene X\*9/\*99 and Y\*77:01"
    assert escape_star_alleles("`*9/*99` kept, *9 escaped") == r"`*9/*99` kept, \*9 escaped"
    assert escape_star_alleles("**1. Heading** and *italic*") == "**1. Heading** and *italic*"
    assert escape_star_alleles("odd ` backtick *9") == "odd ` backtick *9"


def test_star_alleles_are_not_italicised():
    body = markdown_to_html_body("Negative for Z*90:01, Z*91:01, Z*92:01\n")
    assert "<em>" not in body
    assert "Z*90:01, Z*91:01, Z*92:01" in body


def test_orientation_follows_widest_table():
    assert max_table_columns("no tables here") == 0
    assert max_table_columns(_table(3) + "\n\n" + _table(9)) == 9
    assert page_orientation(_table(6)) == "portrait"
    assert page_orientation(_table(7)) == "landscape"


def test_build_html_document_escapes_title_and_sets_page_size():
    doc = build_html_document("a<b", "# Title\n\n" + _table(8), footer="Rendered from x/y.md")
    assert "<title>a&lt;b</title>" in doc
    assert "size: letter landscape" in doc
    assert "<h1>Title</h1>" in doc
    assert 'class="footer">Rendered from x/y.md<' in doc
    assert "size: letter portrait" in build_html_document("t", "plain")


def test_report_markdown_files_lists_top_level_md_only(tmp_path):
    (tmp_path / "b_report.md").write_text("x", encoding="utf-8")
    (tmp_path / "a_card.md").write_text("x", encoding="utf-8")
    (tmp_path / "data.tsv").write_text("x", encoding="utf-8")
    (tmp_path / "pharmcat").mkdir()
    (tmp_path / "pharmcat" / "nested.md").write_text("x", encoding="utf-8")
    assert [p.name for p in report_markdown_files(tmp_path)] == ["a_card.md", "b_report.md"]


def test_find_browser_prefers_env_override():
    assert find_browser({"CHROME": "/opt/custom/chrome"}) == "/opt/custom/chrome"


def test_chrome_print_command_for_native_browser(tmp_path):
    cmd = chrome_print_command("/usr/bin/chromium", tmp_path / "r.html", tmp_path / "r.pdf", tmp_path / "prof")
    assert cmd[0] == "/usr/bin/chromium"
    assert "--headless" in cmd and "--no-pdf-header-footer" in cmd
    assert f"--print-to-pdf={tmp_path / 'r.pdf'}" in cmd
    assert cmd[-1] == (tmp_path / "r.html").resolve().as_uri()


FAKE_BROWSER = """#!{python}
import sys
args = sys.argv[1:]
out = next(a.split("=", 1)[1] for a in args if a.startswith("--print-to-pdf="))
html_url = args[-1]
assert html_url.startswith("file://"), html_url
if "{mode}" == "ok":
    open(out, "wb").write(b"%PDF-1.4 synthetic")
elif "{mode}" == "junk":
    open(out, "wb").write(b"<html>not a pdf")
sys.exit(0 if "{mode}" == "ok" else 3)
"""


def _fake_browser(tmp_path: Path, mode: str) -> str:
    exe = tmp_path / f"fakebrowser_{mode}"
    exe.write_text(FAKE_BROWSER.format(python=sys.executable, mode=mode), encoding="utf-8")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    return str(exe)


@pytest.mark.skipif(sys.platform == "win32", reason="fake browser is a shebang script")
def test_render_pdf_moves_valid_pdf_into_place(tmp_path):
    md = tmp_path / "card.md"
    md.write_text("# Card\n\n" + _table(2), encoding="utf-8")
    pdf = tmp_path / "card.pdf"
    render_pdf(md, pdf, _fake_browser(tmp_path, "ok"), tmp_path / "work")
    assert pdf.read_bytes().startswith(b"%PDF-")
    assert list((tmp_path / "work").iterdir()) == []  # temp HTML and profile cleaned up


@pytest.mark.skipif(sys.platform == "win32", reason="fake browser is a shebang script")
@pytest.mark.parametrize("mode", ["junk", "none"])
def test_render_pdf_keeps_old_pdf_when_browser_fails(tmp_path, mode):
    md = tmp_path / "card.md"
    md.write_text("# Card\n", encoding="utf-8")
    pdf = tmp_path / "card.pdf"
    pdf.write_bytes(b"%PDF-old")
    with pytest.raises(RuntimeError, match="did not write a PDF"):
        render_pdf(md, pdf, _fake_browser(tmp_path, mode), tmp_path / "work")
    assert pdf.read_bytes() == b"%PDF-old"
