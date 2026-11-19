#!/usr/bin/env python3
"""Check every cited reference against CrossRef.

    python3 tools/check_refs.py [manuscript/reference.bib] [main.tex si.tex ...]

Collects the keys cited with \\cite{...}, looks each one up by DOI where the
entry has one and by title otherwise, and reports any disagreement in year,
volume, first page, journal or title.  Abbreviated journal names (ACS style)
are not flagged.  A lookup CrossRef answers with an HTTP error, such as an
unknown DOI, is reported as ERR.  Needs only the standard library and network
access.
"""
import html
import json
import re
import sys
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request

if len(sys.argv) > 1:
    BIB = sys.argv[1]
else:
    BIB = 'manuscript/reference.bib'
if len(sys.argv) > 2:
    TEX = sys.argv[2:]
else:
    TEX = ['manuscript/main.tex', 'manuscript/si.tex']
USER_AGENT = 'ref-check (mailto:qzhang234@anl.gov)'


def field(body, name):
    """Value of one field of a BibTeX entry body, without its braces or quotes;
    '' if the entry has no such field."""
    m = re.search(rf'(?<![a-zA-Z]){name}\s*=\s*', body, re.I)
    if not m:
        return ''
    i = m.end()
    if body[i] == '{':
        depth = 0
        for j in range(i, len(body)):
            if body[j] == '{':
                depth += 1
            elif body[j] == '}':
                depth -= 1
                if depth == 0:
                    return body[i + 1:j]
        return ''
    if body[i] == '"':
        # a quoted value may contain braces; only a quote at brace depth zero closes it
        depth = 0
        for j in range(i + 1, len(body)):
            if body[j] == '{':
                depth += 1
            elif body[j] == '}':
                depth -= 1
            elif body[j] == '"' and depth == 0:
                return body[i + 1:j]
        return ''
    return re.split(r',\s*\n|\n', body[i:])[0].strip().rstrip(',')


def clean(s):
    """Plain ASCII text: HTML tags, LaTeX commands and braces removed."""
    s = html.unescape(re.sub(r'<[^>]+>', '', s))
    s = re.sub(r'\\[a-zA-Z]+\s*', '', s)
    s = re.sub(r'[{}$\\~"\'`^]', '', s)
    s = unicodedata.normalize('NFKD', s).encode('ascii', 'ignore').decode()
    return re.sub(r'\s+', ' ', s).strip()


def norm(s):
    """Lower-case letters and digits only, for comparing titles."""
    return re.sub(r'[^a-z0-9]', '', clean(s).lower())


def get(url, tries=4):
    """JSON from CrossRef.  CrossRef answers bursts with HTTP 429, so those are
    retried after a pause; any other HTTP error is raised."""
    request = urllib.request.Request(url, headers={'User-Agent': USER_AGENT})
    for n in range(tries):
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                return json.load(response)
        except urllib.error.HTTPError as e:
            if e.code != 429 or n == tries - 1:
                raise
            time.sleep(2 ** n)


def first_page(pages):
    return pages.split('-')[0]


with open(BIB) as fh:
    bib = fh.read()
cited = []
for tex_file in TEX:
    with open(tex_file) as fh:
        tex = fh.read()
    for m in re.finditer(r'\\cite\{([^}]*)\}', tex):
        for key in m.group(1).split(','):
            key = key.strip()
            if key and key not in cited:
                cited.append(key)

entries = {}
for m in re.finditer(r'@(\w+)\{([^,]+),(.*?)\n\}', bib, re.S):
    entries[m.group(2).strip()] = (m.group(1).lower(), m.group(3))

print(f'checking {len(cited)} cited references against CrossRef\n')
bad = []
for key in cited:
    entry_type, body = entries[key]
    title = clean(field(body, 'title'))
    if entry_type == 'misc' or not title:
        print(f'SKIP  {key}  (non-article: {entry_type})')
        continue
    doi = clean(field(body, 'doi'))
    try:
        if doi:
            msg = get('https://api.crossref.org/works/' + urllib.parse.quote(doi))['message']
        else:
            query = urllib.parse.urlencode({'query.bibliographic': title, 'rows': 6})
            items = get('https://api.crossref.org/works?' + query)['message']['items']
            msg = None
            for item in items:
                if item.get('type') == 'component':
                    continue
                if norm(item.get('title', [''])[0])[:60] == norm(title)[:60]:
                    msg = item
                    break
            if msg is None:
                print(f'??    {key}  no confident CrossRef match for title')
                bad.append((key, 'no match'))
                continue
    except urllib.error.HTTPError as e:
        print(f'ERR   {key}  {e}')
        bad.append((key, f'lookup failed: {e}'))
        continue

    date_parts = msg.get('issued', {}).get('date-parts') or [[None]]
    cr_year = date_parts[0][0]
    cr_volume = str(msg.get('volume', '') or '')
    # Journals that number articles leave `page` empty and give `article-number`.
    cr_page = str(msg.get('page', '') or msg.get('article-number', '') or '')
    cr_journal = (msg.get('container-title') or [''])[0]
    cr_title = (msg.get('title') or [''])[0]
    bib_year = clean(field(body, 'year'))
    bib_volume = clean(field(body, 'volume'))
    bib_pages = clean(field(body, 'pages')).replace('--', '-')
    bib_journal = clean(field(body, 'journal'))

    problems = []
    if bib_year and cr_year and str(cr_year) != bib_year:
        problems.append(f'year bib={bib_year} crossref={cr_year}')
    if bib_volume and cr_volume and bib_volume != cr_volume:
        problems.append(f'volume bib={bib_volume} crossref={cr_volume}')
    if bib_pages and cr_page:
        if first_page(bib_pages).lstrip('0') != first_page(cr_page).lstrip('0'):
            problems.append(f'pages bib={bib_pages} crossref={cr_page}')
    if norm(title)[:50] != norm(cr_title)[:50]:
        problems.append(f'title differs: crossref="{cr_title[:70]}"')
    jb = norm(bib_journal)
    jc = norm(cr_journal)
    abbreviated = '.' in bib_journal          # "Annu. Rev. Biomed. Eng." is ACS style, not an error
    if jb and jc and not abbreviated and not (jb in jc or jc in jb or len(jb) < 12):
        problems.append(f'journal bib="{bib_journal}" crossref="{cr_journal}"')
    if problems:
        print(f'FLAG  {key}')
        for p in problems:
            print(f'         - {p}')
        bad.append((key, problems))
    else:
        print(f'ok    {key}  {cr_journal} {cr_year}, {cr_volume}, {cr_page}')
    time.sleep(0.4)
print(f'\n{len(bad)} reference(s) flagged')
