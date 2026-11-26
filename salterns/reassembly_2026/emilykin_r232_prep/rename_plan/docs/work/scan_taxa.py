"""Scan text documents (docx/md/pdf/pptx/html/tsv/csv segments) for R220 taxon names that change in R232, R220 species
strings, and GTDB version/release statements. Writes taxa_occurrences_raw.tsv."""
import json, csv, re, sys
from collections import Counter
NS=list(csv.DictReader(open('r220_name_stats.tsv'),delimiter='\t'))
SP=[r for r in csv.DictReader(open('rank_name_changes.tsv'),delimiter='\t') if r['rank']=='s' and r['r220']!='s__']
FILES=set(l.strip() for l in open(sys.argv[1]) if l.strip())
terms=[]
for r in NS:
    if r['rank']=='s': continue
    n=r['r220_name']
    # skip names too generic to scan without lineage prefix? keep all, flag partial
    terms.append(dict(term=n, kind='taxon_'+r['rank'], repl=r['r232_targets'], partial=int(r['n_changed'])<int(r['n_mags_r220']),
                      rx=re.compile(r'(?<![A-Za-z0-9])(?:'+r['rank']+'__)?'+re.escape(n)+r'(?![A-Za-z0-9])')))
for r in SP:
    n=r['r220'][3:]
    terms.append(dict(term=n, kind='species', repl=(r['r232'][3:] or '(no species in R232)')+' ['+r['mags']+']', partial=False,
                      rx=re.compile(r'(?<![A-Za-z0-9])(?:s__)?'+re.escape(n).replace(r'\ ', r'\s+')+r'(?![A-Za-z0-9])')))
    # bare epithet code (e.g. sp000585075 / phosphatis_G) also
    ep=n.split(' ',1)[1]
    if re.match(r'sp\d+$',ep) or '_' in ep:
        terms.append(dict(term=ep, kind='species_epithet', repl=r['r232'][3:]+' ['+r['mags']+']', partial=False,
                          rx=re.compile(r'(?<![A-Za-z0-9])'+re.escape(ep)+r'(?![A-Za-z0-9])')))
VER=re.compile(r'GTDB[- ]?Tk[^.;,)]{0,40}|\bGTDB\b[^.;)]{0,60}|\b[Rr](?:elease)?\s?2(?:02|07|14|20|26|32)\b|release ?2(?:14|20|26|32)|\bv?2\.7\.2\b')
ASVCUE=re.compile(r'\bASVs?\b|amplicon|MiDAS|SILVA|16S rRNA gene amplicon|iterativeIDs?\b(?!.*MAG)',re.I)
GTDBCUE=re.compile(r'GTDB|[pcofgs]__|\bMAGs?\b|CAN_\d_bin|coasm_bin|genome',re.I)
out=[]
for line in open('segments.jsonl'):
    s=json.loads(line)
    if s['file'] not in FILES: continue
    t=s['text']
    if not t: continue
    ctxt=t+' '+(s.get('row') or '')
    for d in terms:
        for mo in d['rx'].finditer(t):
            g=bool(GTDBCUE.search(ctxt)); a=bool(ASVCUE.search(ctxt))
            ctx='GTDB/MAG' if g and not a else 'ASV/16S' if a and not g else 'mixed' if a and g else 'unclear'
            if mo.group(0)[1:3]=='__': ctx='GTDB/MAG (lineage string)'
            out.append(dict(file=s['file'],loc=json.dumps(s['loc'],ensure_ascii=False),state=s['state'],kind=d['kind'],term=d['term'],
                match=mo.group(0),proposed=d['repl'],partial_rename=d['partial'],auto_context=ctx,
                quote=t[max(0,mo.start()-150):mo.end()+150].replace('\t',' ').replace('\n',' ')))
    for mo in VER.finditer(t):
        out.append(dict(file=s['file'],loc=json.dumps(s['loc'],ensure_ascii=False),state=s['state'],kind='version',term='GTDB/version',
            match=mo.group(0),proposed='',partial_rename=False,auto_context='version',
            quote=t[max(0,mo.start()-150):mo.end()+200].replace('\t',' ').replace('\n',' ')))
with open('taxa_occurrences_raw.tsv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(out[0]),delimiter='\t'); w.writeheader(); w.writerows(out)
c=Counter((o['file'],o['kind']) for o in out)
for k,v in sorted(c.items()): print(v,*k)
