#!/usr/bin/env python3
"""Confirm the PNAS significance statement has not drifted between its two homes.

It lives twice: as a section of manuscript/main.tex, and as the standalone page
manuscript/significance.tex that the submission portal wants as its own field.
Two copies drift, so compare them word for word and check the 120-word limit.

    python3 tools/check_significance.py
"""
import os
import re
import sys

MS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'manuscript')


def words(t):
    return [w for w in t.split() if re.search(r'[A-Za-z0-9]', w)]


def grab(path, pattern):
    m = re.search(pattern, open(os.path.join(MS, path)).read(), re.S)
    if not m:
        sys.exit(f'FAIL: could not find the statement in {path}')
    return ' '.join(m.group(1).split())


main = grab('main.tex', r'section\*\{Significance Statement\}(.*?)\n\n')
sig = grab('significance.tex', r'Computational Biology\}(.*?)\\vspace\{1\.4em\}')
# the standalone page carries the title block before the text; drop it
sig = sig.split('}', 1)[-1].strip() if sig.startswith('\\vspace') is False else sig
sig = re.sub(r'^.*?\\vspace\{1\.2em\}\s*', '', sig).strip()

n = len(words(main))
print(f'main.tex                : {n} words (limit 120)')
print(f'significance.tex        : {len(words(sig))} words')
ok = True
if n > 120:
    print('FAIL: over the 120-word limit'); ok = False
if main != sig:
    print('FAIL: the two copies have drifted')
    print(f'  main.tex        : {main[:100]}...')
    print(f'  significance.tex: {sig[:100]}...')
    ok = False
else:
    print('identical in both files')
sys.exit(0 if ok else 1)
