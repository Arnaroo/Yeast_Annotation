#!/usr/bin/env python3
"""
compose_figure1.py
=====================
Stacks the wide panel A (plot_figure1a_schematic.py) above the panel B
block (plot_figure1b_and_suppl_fig3.py) into a single page. Both inputs
are vector PDFs at the same width, so the composite is a pure
translation with no scaling and no rasterisation: nothing in either
panel is resampled, and text in both stays selectable.

The width is checked rather than assumed: if the two inputs ever stop
matching, this script fails instead of silently centring one of them.

Input:  figure1a.pdf, figure1b.pdf
Output: figure1.pdf

Run:
  python3 compose_figure1.py --panel-a figure1a.pdf --panel-b figure1b.pdf --out figure1.pdf
"""

import argparse
import sys

from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import RectangleObject

MAX_WIDTH_PT = 7.5 * 72.0  # journal single-page text width, leaves a margin
GAP_PT = 6.0  # vertical breathing space between the two blocks


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-a", required=True)
    ap.add_argument("--panel-b", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    page_a = PdfReader(args.panel_a).pages[0]
    page_b = PdfReader(args.panel_b).pages[0]
    wa, ha = float(page_a.mediabox.width), float(page_a.mediabox.height)
    wb, hb = float(page_b.mediabox.width), float(page_b.mediabox.height)
    print(f"panel A  {wa:.1f} x {ha:.1f} pt   ({wa / 72:.2f} x {ha / 72:.2f} in)")
    print(f"panel B  {wb:.1f} x {hb:.1f} pt   ({wb / 72:.2f} x {hb / 72:.2f} in)")

    if abs(wa - wb) > 0.5:
        sys.exit(f"FAIL: panel widths differ, {wa:.1f} vs {wb:.1f} pt. Redraw one of them rather than scaling here.")

    width, height = wa, ha + GAP_PT + hb
    print(f"composite {width:.1f} x {height:.1f} pt ({width / 72:.2f} x {height / 72:.2f} in)")
    if width > MAX_WIDTH_PT:
        sys.exit(f"FAIL: composite is {width / 72:.2f} in wide, over the {MAX_WIDTH_PT / 72:.2f} in limit.")

    writer = PdfWriter()
    page = writer.add_blank_page(width=width, height=height)
    # PDF origin is bottom-left. Panel A goes on top, panel B underneath.
    page.merge_transformed_page(page_a, Transformation().translate(0, hb + GAP_PT))
    page.merge_transformed_page(page_b, Transformation().translate(0, 0))
    page.mediabox = RectangleObject((0, 0, width, height))

    with open(args.out, "wb") as fh:
        writer.write(fh)
    print(f"\nwrote {args.out}")
    print("No scaling applied. Both blocks are placed by translation only.")


if __name__ == "__main__":
    main()
