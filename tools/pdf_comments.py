#!/usr/bin/env python3
"""Bind PDF review comments to the text they highlight.

    python3 tools/pdf_comments.py review/comments.xfdf manuscript/build/main.pdf

Reading a commented PDF by eye does not work when the comments are passed on
as screenshots: the highlights and the notes arrive separately and have to be
matched up by hand, one at a time.  This does the matching.

Export the comments from Acrobat first -- Comments pane, the "..." menu,
"Export All to Data File", save as .xfdf.  That file is XML and holds every
note's text plus the page and coordinates of the span it marks.  The words
underneath those coordinates come from `pdftotext -bbox-layout`, which reports
a box for every word on the page.  Intersect the two and each note is printed
with the text it is attached to, in document order.

Only the standard library and poppler's pdftotext are needed; no PDF package
is installed on the beamline machines.
"""
import html
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

# XFDF puts every element in this namespace; strip it rather than write it out.
NS = re.compile(r'\{[^}]*\}')

# A word counts as highlighted when this much of it lies inside the mark.  A
# highlight drawn by hand overshoots the line and clips the first and last
# glyphs, so anything near 1.0 drops the words at both ends of every span.
OVERLAP = 0.45


def annots(xfdf):
    """(page, ymin, quads, note) per markup annotation, in reading order.

    XFDF coordinates are PDF user space: origin bottom left, y increasing
    upwards.  `coords` is a flat list of quadrilaterals, x1,y1...x4,y4 per
    highlighted line, so one annotation spanning two lines has eight pairs.
    `rect` is the fallback for note types that carry no quads.
    """
    out = []
    for el in ET.parse(xfdf).getroot().iter():
        tag = NS.sub('', el.tag)
        if tag not in ('highlight', 'text', 'underline', 'squiggly', 'strikeout'):
            continue
        page = int(el.get('page', 0))
        note = ''
        for child in el:
            if NS.sub('', child.tag) == 'contents':
                note = ''.join(child.itertext()).strip()
        if not note:
            note = html.unescape(el.get('contents', '') or '').strip()

        nums = [float(v) for v in re.split(r'[,\s]+', el.get('coords', '').strip()) if v]
        quads = [nums[i:i + 8] for i in range(0, len(nums) - 7, 8)]
        if not quads:
            r = [float(v) for v in re.split(r'[,\s]+', el.get('rect', '').strip()) if v]
            if len(r) == 4:
                x0, y0, x1, y1 = r
                quads = [[x0, y1, x1, y1, x0, y0, x1, y0]]
        if not quads:
            continue
        boxes = [(min(q[0::2]), min(q[1::2]), max(q[0::2]), max(q[1::2])) for q in quads]
        out.append((page, -max(b[3] for b in boxes), boxes, note))
    return sorted(out, key=lambda a: (a[0], a[1]))


def words(pdf):
    """{page: (height, [(x0, y0, x1, y1, word)])}, y measured DOWN from the top."""
    xml = subprocess.run(['pdftotext', '-bbox-layout', pdf, '-'],
                         capture_output=True, text=True, check=True).stdout
    pages, page = {}, -1
    for m in re.finditer(r'<page width="([\d.]+)" height="([\d.]+)"'
                         r'|<word xMin="([\d.]+)" yMin="([\d.]+)" '
                         r'xMax="([\d.]+)" yMax="([\d.]+)">(.*?)</word>', xml, re.S):
        if m.group(1):
            page += 1
            pages[page] = (float(m.group(2)), [])
        elif page >= 0:
            pages[page][1].append((float(m.group(3)), float(m.group(4)),
                                   float(m.group(5)), float(m.group(6)),
                                   html.unescape(m.group(7))))
    return pages


def covered(boxes, height, page_words):
    """Words whose area overlaps any quad by more than OVERLAP."""
    hit = []
    for x0, ytop, x1, ybot, w in page_words:
        # flip pdftotext's top-down y into the bottom-up y that XFDF uses
        wy0, wy1 = height - ybot, height - ytop
        area = max(x1 - x0, 1e-6) * max(wy1 - wy0, 1e-6)
        for bx0, by0, bx1, by1 in boxes:
            ox = min(x1, bx1) - max(x0, bx0)
            oy = min(wy1, by1) - max(wy0, by0)
            if ox > 0 and oy > 0 and ox * oy / area >= OVERLAP:
                hit.append(w)
                break
    return ' '.join(hit)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    xfdf, pdf = sys.argv[1], sys.argv[2]
    pages = words(pdf)
    for n, (page, _, boxes, note) in enumerate(annots(xfdf), 1):
        height, page_words = pages.get(page, (0, []))
        text = covered(boxes, height, page_words) or '(no text under this mark)'
        print(f'--- C{n}  page {page + 1}')
        print(f'    marked : {text}')
        print(f'    comment: {note or "(highlight only, no note)"}\n')


if __name__ == '__main__':
    main()
