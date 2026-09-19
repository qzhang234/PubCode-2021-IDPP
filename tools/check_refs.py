#!/usr/bin/env python3
"""Check every cited reference against CrossRef.

    python3 tools/check_refs.py [manuscript/reference.bib] [main.tex si.tex ...]

Collects the keys actually cited, looks each one up by DOI where the entry
has one and by title otherwise, and reports any disagreement in year,
volume, first page, journal or title.  Journal abbreviations are not
flagged: ACS style uses them.  Needs only the standard library.
"""
import html, json, re, sys, time, urllib.parse, urllib.request, unicodedata

BIB=sys.argv[1] if len(sys.argv)>1 else 'manuscript/reference.bib'
TEX=sys.argv[2:] or ['manuscript/main.tex','manuscript/si.tex']
UA='ref-check (mailto:qzhang234@anl.gov)'

def field(body, name):
    m = re.search(rf'(?<![a-zA-Z]){name}\s*=\s*', body, re.I)
    if not m: return ''
    i = m.end()
    if body[i] == '{':
        d=0; j=i
        while j < len(body):
            if body[j]=='{': d+=1
            elif body[j]=='}':
                d-=1
                if d==0: return body[i+1:j]
            j+=1
        return ''
    if body[i] == '"':
        # BibTeX allows braces inside a quoted value, so only a quote at brace
        # depth zero closes it.  Treating '"' as its own opener (as before) made
        # the closing quote look like a second opener and returned nothing.
        d=0; j=i+1
        while j < len(body):
            if body[j]=='{': d+=1
            elif body[j]=='}': d-=1
            elif body[j]=='"' and d==0: return body[i+1:j]
            j+=1
        return ''
    return re.split(r',\s*\n|\n', body[i:])[0].strip().rstrip(',')

def clean(s):
    s = html.unescape(re.sub(r'<[^>]+>', '', s))
    s = re.sub(r'\\[a-zA-Z]+\s*', '', s)
    s = re.sub(r'[{}$\\~"\'`^]', '', s)
    s = unicodedata.normalize('NFKD', s).encode('ascii','ignore').decode()
    return re.sub(r'\s+', ' ', s).strip()

def norm(s):
    return re.sub(r'[^a-z0-9]', '', clean(s).lower())

def get(url, tries=4):
    """CrossRef throttles bursts with 429; back off rather than lose the entry."""
    req = urllib.request.Request(url, headers={'User-Agent': UA})
    for n in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code != 429 or n == tries - 1:
                raise
            time.sleep(2 ** n)

src = open(BIB).read()
cited=[]
for t in [open(f).read() for f in TEX]:
    for m in re.finditer(r'\\cite\{([^}]*)\}', t):
        for k in m.group(1).split(','):
            k=k.strip()
            if k and k not in cited: cited.append(k)

ents={}
for m in re.finditer(r'@(\w+)\{([^,]+),(.*?)\n\}', src, re.S):
    ents[m.group(2).strip()]=(m.group(1).lower(), m.group(3))

print(f'checking {len(cited)} cited references against CrossRef\n')
bad=[]
for k in cited:
    typ, body = ents[k]
    title = clean(field(body,'title'))
    if typ == 'misc' or not title:
        print(f'SKIP  {k}  (non-article: {typ})'); continue
    doi = clean(field(body,'doi'))
    try:
        if doi:
            msg = get('https://api.crossref.org/works/'+urllib.parse.quote(doi))['message']
        else:
            q = urllib.parse.urlencode({'query.bibliographic': title, 'rows': 6})
            items = get('https://api.crossref.org/works?'+q)['message']['items']
            items = [i for i in items if i.get('type') != 'component']
            msg = next((i for i in items if norm(i.get('title',[''])[0])[:60]==norm(title)[:60]), None)
            if msg is None:
                print(f'??    {k}  no confident CrossRef match for title'); bad.append((k,'no match')); continue
    except Exception as e:
        print(f'ERR   {k}  {e}'); bad.append((k,f'lookup failed: {e}')); continue

    cy = (msg.get('issued',{}).get('date-parts') or [[None]])[0][0]
    cv = str(msg.get('volume','') or '')
    # Journals that number articles rather than paginate leave `page` empty and
    # put the number in `article-number`.  Without this the page check silently
    # passed for a third of the bibliography -- it is how a "1--10" placeholder
    # survived in an entry whose real locator is article 196.
    cp = str(msg.get('page','') or msg.get('article-number','') or '')
    cj = (msg.get('container-title') or [''])[0]
    ct = (msg.get('title') or [''])[0]
    by, bv = clean(field(body,'year')), clean(field(body,'volume'))
    bp, bj = clean(field(body,'pages')).replace('--','-'), clean(field(body,'journal'))

    probs=[]
    if by and cy and str(cy)!=by: probs.append(f'year bib={by} crossref={cy}')
    if bv and cv and bv!=cv:      probs.append(f'volume bib={bv} crossref={cv}')
    if bp and cp and bp.split('-')[0]!=cp.split('-')[0].lstrip('0'):
        if bp.split('-')[0].lstrip('0')!=cp.split('-')[0].lstrip('0'):
            probs.append(f'pages bib={bp} crossref={cp}')
    if norm(title)[:50]!=norm(ct)[:50]: probs.append(f'title differs: crossref="{ct[:70]}"')
    jb, jc = norm(bj), norm(cj)
    abbrev = '.' in bj            # "Annu. Rev. Biomed. Eng." is ACS style, not an error
    if jb and jc and not abbrev and not (jb in jc or jc in jb or len(jb) < 12):
        probs.append(f'journal bib="{bj}" crossref="{cj}"')
    if probs:
        print(f'FLAG  {k}'); [print(f'         - {p}') for p in probs]; bad.append((k,probs))
    else:
        print(f'ok    {k}  {cj} {cy}, {cv}, {cp}')
    time.sleep(0.4)
print(f'\n{len(bad)} reference(s) flagged')
