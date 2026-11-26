import zipfile,re,json,sys
from lxml import etree
W='http://schemas.openxmlformats.org/wordprocessingml/2006/main'
NS={'w':W}
def accepted_paras(path, part='word/document.xml'):
    z=zipfile.ZipFile(path); x=etree.fromstring(z.read(part))
    for d in x.xpath('//w:del|//w:moveFrom',namespaces=NS):
        d.getparent().remove(d)
    return x
def para_text(p): return ''.join(t.text or '' for t in p.xpath('.//w:t',namespaces=NS))
def tables(x):
    out=[]
    for t in x.xpath('//w:tbl',namespaces=NS):
        rows=[]
        for tr in t.xpath('./w:tr',namespaces=NS):
            rows.append([' '.join(para_text(p) for p in tc.xpath('./w:p',namespaces=NS)) for tc in tr.xpath('./w:tc',namespaces=NS)])
        out.append(rows)
    return out
if __name__=='__main__':
    m=json.load(open('/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan/apply_scratch/verify_ids/ref_map.json'));NEW=m['new'];OLD=m['old']
    sets=json.load(open('/home/freiburger/Documents/MicrobiomeNotebooks/salterns/reassembly_2026/emilykin_r232_prep/rename_plan/apply_scratch/verify_ids/id_sets.json'))
    binrx=re.compile(r'(CAN_[1-5]_bin\.\d+|coasm_bin\.\d+)')
    for f in sys.argv[1:]:
        x=accepted_paras(f)
        print('==',f)
        for ti,t in enumerate(tables(x)):
            nb=0;bad=[]
            for ri,row in enumerate(t):
                bins=[c for c in row if binrx.fullmatch(c.strip())]
                if len(bins)!=1: continue
                b=bins[0].strip(); nb+=1
                if b not in NEW: bad.append((ri,b,'not291')); continue
                ids=[c.strip() for c in row if re.fullmatch(r'[A-Za-z0-9_\-]+\.\d+(_m)?',c.strip()) and not binrx.fullmatch(c.strip())]
                if ids and NEW[b] not in ids: bad.append((ri,b,ids,NEW[b]))
            if nb: print(' table',ti,'rows with bin',nb,'bad',bad)
        text='\n'.join(para_text(p) for p in x.xpath('//w:p',namespaces=NS))
        hits=[i for i in sets['old_only_mag']+sets['old_reused'] if re.search(r'(?<![A-Za-z0-9])'+re.escape(i)+r'(?![0-9A-Za-z])',text)]
        print(' old MAG-only _m IDs in accepted text:',hits)
        inh=[(i,len(re.findall(r'(?<![A-Za-z0-9])'+re.escape(i)+r'(?![0-9A-Za-z])',text))) for i in sets['old_inh']]
        print(' inherited old IDs (count) in accepted text:',[x for x in inh if x[1]])
        # ID near bin mismatch
        for mm in re.finditer(r'([A-Za-z][A-Za-z0-9_\-]*\.\d+(?:_m)?)[^\n]{0,6}?(CAN_[1-5]_bin\.\d+|coasm_bin\.\d+)',text):
            i,b=mm.group(1),mm.group(2)
            if b in NEW and i!=NEW[b] and (i in OLD.values() or i in NEW.values()):
                print('  near-bin mismatch:',i,b,'R232:',NEW[b],'|',text[max(0,mm.start()-60):mm.end()+40].replace('\n',' '))
        for mm in re.finditer(r'(CAN_[1-5]_bin\.\d+|coasm_bin\.\d+)[^\n]{0,6}?([A-Za-z][A-Za-z0-9_\-]*\.\d+(?:_m)?)(?![0-9A-Za-z])',text):
            b,i=mm.group(1),mm.group(2)
            if b in NEW and i!=NEW[b] and (i in OLD.values() or i in NEW.values()):
                print('  near-bin mismatch:',b,i,'R232:',NEW[b],'|',text[max(0,mm.start()-60):mm.end()+40].replace('\n',' '))
