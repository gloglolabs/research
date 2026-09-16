"""Build the two-page PDF from research-note.md using ReportLab."""

from html import escape
from pathlib import Path
import re

from reportlab import rl_config
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

ROOT = Path(__file__).resolve().parent
INK = colors.HexColor("#152A38")
MUTED = colors.HexColor("#53616C")


def inline(text):
    text = text.replace(r"\*", "*")
    parts = re.split(r"(\*\*.+?\*\*|`[^`]+`|\[[^\]]+\]\([^)]+\))", text)
    output = []
    for part in parts:
        if part.startswith("**") and part.endswith("**"):
            output.append(f"<b>{escape(part[2:-2])}</b>")
        elif part.startswith("`") and part.endswith("`"):
            output.append(f'<font name="Courier" size="9.1">{escape(part[1:-1])}</font>')
        elif part.startswith("[") and "](" in part:
            label, url = part[1:].split("](", 1)
            output.append(f'<link href="{escape(url[:-1], quote=True)}" color="#1B566D"><u>{escape(label)}</u></link>')
        else:
            output.append(escape(part))
    return "".join(output)


def footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#D6DEE3"))
    canvas.line(44, 38, letter[0] - 44, 38)
    canvas.setFont("NoteSans", 8)
    canvas.setFillColor(MUTED)
    canvas.drawString(44, 25, "GloGlo Labs | CHGA-INS research note | v1.1")
    canvas.drawRightString(letter[0] - 44, 25, str(doc.page))
    canvas.restoreState()


def main():
    font_dir = next(Path(p) for p in rl_config.TTFSearchPath if (Path(p) / "Vera.ttf").is_file())
    pdfmetrics.registerFont(TTFont("NoteSans", str(font_dir / "Vera.ttf")))
    pdfmetrics.registerFont(TTFont("NoteSans-Bold", str(font_dir / "VeraBd.ttf")))
    pdfmetrics.registerFontFamily("NoteSans", normal="NoteSans", bold="NoteSans-Bold")
    body = ParagraphStyle("Body", fontName="NoteSans", fontSize=9.7, leading=13.3, textColor=INK, spaceAfter=8)
    title = ParagraphStyle("Title", parent=body, fontName="NoteSans-Bold", fontSize=18, leading=22, spaceAfter=8)
    heading = ParagraphStyle("Heading", parent=body, fontName="NoteSans-Bold", fontSize=11, leading=14, spaceBefore=8, spaceAfter=6, keepWithNext=True)
    meta = ParagraphStyle("Meta", parent=body, fontSize=7.7, leading=11, textColor=MUTED, spaceAfter=12)
    cell = ParagraphStyle("Cell", parent=body, fontSize=8.5, leading=11, spaceAfter=0)
    number = ParagraphStyle("Number", parent=cell, alignment=TA_CENTER)
    doc = SimpleDocTemplate(str(ROOT / "chga-ins-research-note.pdf"), pagesize=letter, leftMargin=44, rightMargin=44, topMargin=37, bottomMargin=48, title="L9C and T17C in the CHGA leader of engineered insulin", author="GloGlo Labs", subject="Computational follow-up to Kobaisi et al. (2026)")
    story = []
    for block in (ROOT / "research-note.md").read_text().strip().split("\n\n"):
        if block == "<!-- pagebreak -->":
            story.append(PageBreak())
        elif block.startswith("# "):
            story.append(Paragraph(inline(block[2:]), title))
        elif block.startswith("## "):
            story.append(Paragraph(inline(block[3:]), heading))
        elif block.startswith("|"):
            rows = [[c.strip() for c in line.strip().strip("|").split("|")] for line in block.splitlines()]
            rows = [row for row in rows if not all(re.fullmatch(r":?-+:?", c) for c in row)]
            data = [[Paragraph(("<b>" if i == 0 else "") + inline(value) + ("</b>" if i == 0 else ""), cell if j == 0 else number) for j, value in enumerate(row)] for i, row in enumerate(rows)]
            table = Table(data, colWidths=[186, 65, 65, 65, 143], hAlign="LEFT")
            table.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#EAF0F3")), ("LINEBELOW", (0, 0), (-1, 0), .5, colors.HexColor("#B5C4CD")), ("LINEBELOW", (0, 1), (-1, -1), .3, colors.HexColor("#D9E1E6")), ("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("LEFTPADDING", (0, 0), (-1, -1), 7), ("RIGHTPADDING", (0, 0), (-1, -1), 7), ("TOPPADDING", (0, 0), (-1, -1), 7), ("BOTTOMPADDING", (0, 0), (-1, -1), 7)]))
            story.extend([table, Spacer(1, 9)])
        else:
            text = "<br/>".join(inline(" ".join(part.splitlines())) for part in re.split(r" {2,}\n|\\\n", block))
            story.append(Paragraph(text, meta if block.startswith("GloGlo Labs |") else body))
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    print(ROOT / "chga-ins-research-note.pdf")


if __name__ == "__main__":
    main()
