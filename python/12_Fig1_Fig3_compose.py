#!/usr/bin/env python3
"""
12 - Figures 1 and 3: places the panels drawn by R/10_Fig1_draw.R and R/11_Fig3d_draw.R onto the template pages

assets/Fig1_template.pdf is the Figure 1 page without the data of
panels f and g (cleared completely), the marks of panel h (error bars, points, interaction-test circle) and the
discovery numbers of panel a. assets/Fig3_template.pdf is the Figure 3 page without the marks and the HR labels
of panel d. The panels of scripts 10 and 11 are transparent A4 pages drawn 16 times larger and are scaled onto the
templates.

Usage : python3 python/12_Fig1_Fig3_compose.py
        python3 python/12_Fig1_Fig3_compose.py --make-templates "Fig 1_revised.pdf" "Fig 3.pdf"
                (rebuilds the templates from the figure pages)
Output: ./output/fig1/Fig1.pdf and Fig1.png, ./output/fig3/Fig3.pdf and Fig3.png (300 dpi, figure area)
Needs : Python 3.8+, PyMuPDF
"""
import argparse
import os
import re
import sys

import fitz

HERE = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.join(HERE, "..", "assets")
FIGS = {
    "fig1": dict(template=os.path.join(ASSETS, "Fig1_template.pdf"), panels="./output/fig1/Fig1_panels.pdf",
                 out="./output/fig1/Fig1.pdf", clip=(10, 20, 585, 372)),
    "fig3": dict(template=os.path.join(ASSETS, "Fig3_template.pdf"), panels="./output/fig3/Fig3_panels.pdf",
                 out="./output/fig3/Fig3.pdf", clip=(5, 5, 430, 170)),
}
# Figure 1, panel a: the discovery number strings of the published page (emptied in the content stream, so that
# nothing else moves; script 10 writes the numbers)
FIG1_A_STRINGS = ["n = 687, 189 DDFS events", "n = 415", "103 events", "n = 114", "36 events", "n = 158", "50 events",
                  "n = 360, 77 events", "n = 152", "51 events", "n = 208", "26 events"]
# Figure 1, panels f and g: cleared completely (text and line art, then a white cover), without the panel letters
FIG1_CLEAR = [(187.5, 266.0, 300.0, 362.5), (309.5, 266.0, 422.0, 362.5)]
# Figure 1, panel h, and Figure 3, panel d: only line art lying completely inside these areas (the marks) is removed.
# Axes, grid lines, frames and labels reach outside them and stay.
FIG1_MARKS = [(463.0, 259.0, 512.0, 268.0), (516.5, 259.0, 540.5, 268.0)]
FIG3_MARKS = [(295.0, 34.0, 331.5, 61.5)]
FIG3_TEXT = [(334.5, 33.0, 420.0, 60.5)]                       # panel d, HR labels


def drop_private(doc, page):
    for key in ("PieceInfo", "LastModified", "Thumb"):     # Illustrator private data and the stale thumbnail
        doc.xref_set_key(page.xref, key, "null")


def remove_marks(page, rects):
    for r in rects:
        page.add_redact_annot(fitz.Rect(r), fill=False)
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED,
                          text=fitz.PDF_REDACT_TEXT_NONE)


def make_fig1(published):
    doc = fitz.open(published)
    page = doc[0]
    hits = 0
    for x in page.get_contents():
        raw = doc.xref_stream(x)
        for s in FIG1_A_STRINGS:
            pat = re.escape(b"(" + s.encode("latin1") + b")Tj")
            n = len(re.findall(pat, raw))
            if n:
                raw = re.sub(pat, b"()Tj", raw)
                hits += n
        doc.update_stream(x, raw)
    if hits != len(FIG1_A_STRINGS):
        sys.exit(f"panel a: {hits} of {len(FIG1_A_STRINGS)} number strings found, template not written")
    for r in FIG1_CLEAR:
        page.add_redact_annot(fitz.Rect(r), fill=(1, 1, 1))
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED,
                          text=fitz.PDF_REDACT_TEXT_REMOVE)
    for r in FIG1_CLEAR:                                     # white cover for form XObjects (confidence bands)
        page.draw_rect(fitz.Rect(r), color=None, fill=(1, 1, 1), overlay=True)
    remove_marks(page, FIG1_MARKS)
    drop_private(doc, page)
    left = [s for s in FIG1_A_STRINGS if s in page.get_text("text")]
    if left:
        sys.exit(f"panel a strings still on the page: {left}")
    doc.save(FIGS["fig1"]["template"], garbage=4, deflate=True)
    print("template written:", FIGS["fig1"]["template"])


def make_fig3(published):
    doc = fitz.open(published)
    page = doc[0]
    for r in FIG3_TEXT:
        page.add_redact_annot(fitz.Rect(r), fill=False)
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                          text=fitz.PDF_REDACT_TEXT_REMOVE)
    remove_marks(page, FIG3_MARKS)
    drop_private(doc, page)
    doc.save(FIGS["fig3"]["template"], garbage=4, deflate=True)
    print("template written:", FIGS["fig3"]["template"])


def compose(key):
    f = FIGS[key]
    doc = fitz.open(f["template"])
    page = doc[0]
    panels = fitz.open(f["panels"])
    page.show_pdf_page(page.rect, panels, 0, keep_proportion=False, overlay=True)
    os.makedirs(os.path.dirname(f["out"]), exist_ok=True)
    doc.save(f["out"], garbage=4, deflate=True)
    fitz.open(f["out"])[0].get_pixmap(dpi=300, clip=fitz.Rect(f["clip"])).save(f["out"].replace(".pdf", ".png"))
    print("written:", f["out"])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--make-templates", nargs=2, metavar=("FIG1_PDF", "FIG3_PDF"))
    ap.add_argument("figures", nargs="*", default=["fig1", "fig3"])
    a = ap.parse_args()
    if a.make_templates:
        os.makedirs(ASSETS, exist_ok=True)
        make_fig1(a.make_templates[0])
        make_fig3(a.make_templates[1])
    else:
        for k in a.figures:
            compose(k)
