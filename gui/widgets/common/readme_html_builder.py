#!/usr/bin/env python3
# readme_html_builder.py — Server-side Python Markdown to Interactive HTML converter for README.md.
# PhyloPhere | gui/widgets/common/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

import html
import json
import re
from pathlib import Path


def parse_inline(s: str) -> str:
    """Helper to convert inline markdown (code, bold, italic, links)."""
    # inline code
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    # bold
    s = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", s)
    # italic
    s = re.sub(r"\*([^*]+)\*", r"<em>\1</em>", s)
    # links
    s = re.sub(
        r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2" target="_blank">\1</a>', s
    )
    return s


def parse_table_block(table_text: str) -> str:
    """Helper to convert markdown table into clean HTML table."""
    lines = [l.strip() for l in table_text.strip().split("\n") if l.strip()]
    if len(lines) < 2:
        return f"<pre>{html.escape(table_text)}</pre>"

    headers = [c.strip() for c in lines[0].strip("|").split("|")]
    rows = []
    for line in lines[2:]:
        cols = [c.strip() for c in line.strip("|").split("|")]
        rows.append(cols)

    th_html = "".join(
        [
            f'<th bgcolor="#f6f8fa" style="padding: 6px 10px; border: 1px solid #d0d7de; text-align: left;">{parse_inline(h)}</th>'
            for h in headers
        ]
    )
    tr_html = ""
    for idx, row in enumerate(rows):
        bg = ' bgcolor="#f6f8fa"' if idx % 2 == 1 else ""
        tds = "".join(
            [
                f'<td{bg} style="padding: 6px 10px; border: 1px solid #d0d7de;">{parse_inline(c)}</td>'
                for c in row
            ]
        )
        tr_html += f"<tr>{tds}</tr>\n"

    return f'<div class="table-container" style="overflow-x: auto; margin: 8px 0;"><table border="1" cellspacing="0" cellpadding="4" style="border-collapse: collapse; width: 100%; font-size: 13px;"><thead><tr>{th_html}</tr></thead><tbody>{tr_html}</tbody></table></div>'


def build_qtextbrowser_html(
    md_text: str, expanded_sections: set[int]
) -> tuple[str, int]:
    """Converts README.md into rich HTML specifically formatted for QTextBrowser:
    - Code blocks extracted first so code comments (# ...) are never misparsed as headings
    - Replaces ![logo](src) with explicit <img height="40"> (preserving natural aspect ratio!)
    - Headings (h1 through h6) generated with <a name="slug"> anchors for clickable TOC navigation
    - Replaces <details><summary> with clickable <a href="toggle:SEC_ID"> toggle anchors
    - Returns (rendered_html, total_details_count)
    """

    # 1. PRESERVE CODE BLOCKS FIRST (prevents code comments like # 1. set up... from being turned into <h1>)
    code_blocks = []

    def save_code(m):
        code_txt = html.escape(m.group(2).strip())
        idx = len(code_blocks)
        code_blocks.append(
            f'<pre style="background-color: #f6f8fa; border: 1px solid #d0d7de; border-radius: 4px; padding: 10px; font-family: monospace; font-size: 12px; overflow-x: auto;"><code>{code_txt}</code></pre>'
        )
        return f"___CODE_BLOCK_{idx}___"

    text = re.sub(r"```(\w*)\n([\s\S]*?)```", save_code, md_text)

    # 2. Image Replacement (force small height=40 without width distortion)
    def parse_img(m):
        alt = m.group(1)
        src = m.group(2)
        is_logo = "logo" in alt.lower() or "logo" in src.lower()
        if is_logo or src.endswith("logo.png"):
            return f'<p align="left"><img src="{src}" height="240" alt="{html.escape(alt)}"></p>'
        return f'<p><img src="{src}" alt="{html.escape(alt)}" style="max-width: 100%; height: auto;"></p>'

    text = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)", parse_img, text)

    # 3. Details & Summary replacement
    details_counter = 0

    def replace_details(m):
        nonlocal details_counter
        sec_id = details_counter
        details_counter += 1

        summary_txt = m.group(1).strip()
        body_txt = m.group(2).strip()

        # Parse tables inside details body
        body_html = re.sub(
            r"(\|[^\n]+\|\n\|[-:\s|]+\|\n(?:\|[^\n]+\|\n?)+)",
            lambda tm: parse_table_block(tm.group(1)),
            body_txt,
        )

        # Convert formatting in body
        body_lines = []
        for line in body_html.split("\n"):
            line_str = line.strip()
            if line_str.startswith("- ") or line_str.startswith("* "):
                body_lines.append(f"<li>{parse_inline(line_str[2:])}</li>")
            elif line_str and not line_str.startswith("<"):
                body_lines.append(f"<p>{parse_inline(line_str)}</p>")
            else:
                body_lines.append(line_str)
        body_formatted = "\n".join(body_lines)

        is_open = sec_id in expanded_sections
        arrow = "▼" if is_open else "►"

        if is_open:
            return (
                f'<div style="background-color: #f6f8fa; border: 1px solid #d0d7de; border-radius: 6px; padding: 10px; margin: 12px 0;">'
                f'<a href="toggle:{sec_id}" style="color: #0969da; text-decoration: none;"><b>{arrow} {parse_inline(summary_txt)}</b></a>'
                f'<div style="margin-top: 10px; padding-top: 8px; border-top: 1px solid #d0d7de;">{body_formatted}</div>'
                f"</div>"
            )
        else:
            return (
                f'<div style="background-color: #f6f8fa; border: 1px solid #d0d7de; border-radius: 6px; padding: 10px; margin: 12px 0;">'
                f'<a href="toggle:{sec_id}" style="color: #0969da; text-decoration: none;"><b>{arrow} {parse_inline(summary_txt)}</b></a>'
                f"</div>"
            )

    text = re.sub(
        r"<details>\s*<summary>([\s\S]*?)</summary>([\s\S]*?)</details>",
        replace_details,
        text,
        flags=re.IGNORECASE,
    )

    # 4. Headings (h1 through h6) with anchors for TOC navigation
    slug_counts = {}

    def parse_heading(m):
        hashes, title = m.group(1), m.group(2).strip()
        level = len(hashes)
        slug = re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")
        if slug in slug_counts:
            slug_counts[slug] += 1
            slug = f"{slug}-{slug_counts[slug]}"
        else:
            slug_counts[slug] = 1

        border = (
            ' style="border-bottom: 1px solid #d0d7de; padding-bottom: 4px;"'
            if level <= 2
            else ""
        )
        return f'<a name="{slug}"></a><h{level} id="{slug}"{border}>{parse_inline(title)}</h{level}>'

    text = re.sub(r"^\s*(#{1,6})\s+(.+)$", parse_heading, text, flags=re.MULTILINE)

    # 5. Tables outside details
    text = re.sub(
        r"(\|[^\n]+\|\n\|[-:\s|]+\|\n(?:\|[^\n]+\|\n?)+)",
        lambda tm: parse_table_block(tm.group(1)),
        text,
    )

    # 6. Paragraphs & Lists (do NOT wrap HTML block elements like <h1-h6>, <pre>, <div>, <details>, <table> in <p>)
    out_lines = []
    in_ul = False
    block_tags = (
        "<h1",
        "<h2",
        "<h3",
        "<h4",
        "<h5",
        "<h6",
        "<div",
        "<table",
        "<pre",
        "<details",
        "<a name=",
        "___CODE_BLOCK_",
    )

    for line in text.split("\n"):
        line_str = line.strip()
        if line_str.startswith("- ") or line_str.startswith("* "):
            if not in_ul:
                out_lines.append("<ul>")
                in_ul = True
            out_lines.append(f"<li>{parse_inline(line_str[2:])}</li>")
        else:
            if in_ul:
                out_lines.append("</ul>")
                in_ul = False

            if not line_str:
                continue

            if any(line_str.startswith(tag) for tag in block_tags):
                out_lines.append(line_str)
            else:
                out_lines.append(f"<p>{parse_inline(line_str)}</p>")

    if in_ul:
        out_lines.append("</ul>")

    text = "\n".join(out_lines)

    # 7. Restore Code Blocks LAST
    for idx, cb_html in enumerate(code_blocks):
        text = text.replace(f"___CODE_BLOCK_{idx}___", cb_html)

    full_html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    color: #1f2328;
    background-color: #ffffff;
    line-height: 1.5;
    font-size: 14px;
    padding: 16px;
}}
h1, h2, h3, h4, h5, h6 {{ color: #1f2328; font-weight: 600; }}
a {{ color: #0969da; text-decoration: none; }}
a:hover {{ text-decoration: underline; }}
code {{ background-color: #f6f8fa; padding: 2px 4px; font-family: monospace; font-size: 12px; border-radius: 3px; }}
th, td {{ border: 1px solid #d0d7de; padding: 6px 10px; }}
table {{ border-collapse: collapse; width: 100%; margin: 12px 0; }}
</style>
</head>
<body>
{text}
</body>
</html>"""

    return full_html, details_counter
