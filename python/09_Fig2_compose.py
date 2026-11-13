#!/usr/bin/env python3
"""
09 - Figure 2: places the panels drawn by R/08_Fig2_draw.R onto the template page

assets/Fig2_template.pdf is the Figure 2 page (A4) with the data panels removed, so it holds only the
panel letters, the DISCOVERY COHORT label and the cell schematic of panel a. The panels of script 08 are one
transparent page drawn 16 times larger; it is scaled onto the template.

Usage : python3 python/09_Fig2_compose.py
        python3 python/09_Fig2_compose.py --make-template "Fig 2.pdf"   (rebuilds the template from the figure page)
Output: ./output/fig2/Fig2.pdf and Fig2.png (300 dpi, figure area)
Needs : Python 3.8+, PyMuPDF
"""
import argparse
import os

import fitz

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(HERE, "..", "assets", "Fig2_template.pdf")
PANELS = "./output/fig2/Fig2_panels.pdf"
OUT = "./output/fig2/Fig2.pdf"
# data panel areas of the published page, without the panel letters (x0, y0, x1, y1 in pt from the top left)
CLEAR = [
    (16.5, 70.5, 135.5, 146.5),                                  # a, scatter (the schematic above is kept)
    (144.5, 29.5, 266.5, 70.8),                                  # b
    (138.5, 79.6, 266.5, 146.5),                                 # c
    (275.5, 24.0, 391.5, 146.5), (268.8, 37.0, 275.5, 146.5),    # d
    (398.2, 24.0, 572.0, 146.5), (392.0, 37.0, 398.2, 146.5),    # e
    (27.2, 164.5, 197.5, 292.0), (16.5, 178.3, 27.2, 292.0),     # f
    (204.6, 164.5, 382.5, 292.0), (197.8, 178.3, 204.6, 292.0),  # g
    (389.1, 162.5, 572.0, 292.0), (382.5, 176.7, 389.1, 292.0),  # h
]


def make_template(published):
    doc = fitz.open(published)
    page = doc[0]
    for r in CLEAR:
        page.add_redact_annot(fitz.Rect(r), fill=(1, 1, 1))
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED,
                          text=fitz.PDF_REDACT_TEXT_REMOVE)
    for r in CLEAR:                      # white cover for line art that only partly lies inside an area
        page.draw_rect(fitz.Rect(r), color=None, fill=(1, 1, 1), overlay=True)
    os.makedirs(os.path.dirname(TEMPLATE), exist_ok=True)
    doc.save(TEMPLATE, garbage=4, deflate=True)
    print("template written:", TEMPLATE)


def compose():
    doc = fitz.open(TEMPLATE)
    page = doc[0]
    panels = fitz.open(PANELS)
    page.show_pdf_page(page.rect, panels, 0, keep_proportion=False, overlay=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    doc.save(OUT, garbage=4, deflate=True)
    fitz.open(OUT)[0].get_pixmap(dpi=300, clip=fitz.Rect(10, 15, 520, 296)).save(OUT.replace(".pdf", ".png"))
    print("written:", OUT)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--make-template", metavar="PUBLISHED_FIG2_PDF")
    a = ap.parse_args()
    if a.make_template:
        make_template(a.make_template)
    else:
        compose()
