"""Minimal OOXML tracked-change helper for the GTDB R232 docs stage (surgical, run-level; no paragraph rebuilds).

- Visible ("final view") text of a paragraph = w:t text of runs not inside a w:del.  A match is located in that
  text, the runs at its two ends are split at the match boundaries (rPr cloned), every matched run is wrapped in a
  w:del (w:t -> w:delText) and the replacement is added as new run(s) inside a w:ins right after it.
- A matched run inside an existing w:ins of another author (e.g. the Oct-19 'Claude (A. Freiburger)' Table 1
  insertion) is first isolated: that w:ins is split into before / matched / after parts (each keeps the original
  author/date, new ids); the w:del is nested inside the matched part (schema-legal; Word writes the same) and the new
  w:ins is placed at paragraph level after it, so no text is attributed to the wrong author.
- Comments use the docx skill's comment.add_comment(); range markers are placed at paragraph level around the change.
- Images: the drawing run is wrapped in w:del and a clone pointing at a new media part is inserted in a w:ins."""
import copy, re, sys, os, struct
from lxml import etree

SKILL = "/home/freiburger/.claude/skills/synced/49be370f-b8a1-4df3-8125-e9d4302adaab_67b17de6-c8a0-459f-85ef-8dd13a27989a/docx/scripts"
sys.path.insert(0, SKILL)
from comment import add_comment  # noqa: E402

WNS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
W = "{%s}" % WNS
RNS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
XMLSPACE = "{http://www.w3.org/XML/1998/namespace}space"
AUTHOR = "Claude (GTDB R232 update)"
INITIALS = "CR"


class Doc:
    def __init__(self, unpacked, date):
        self.dir = unpacked
        self.path = os.path.join(unpacked, "word", "document.xml")
        self.tree = etree.parse(self.path)
        self.root = self.tree.getroot()
        self.body = self.root.find(W + "body")
        self.date = date
        ids = [int(e.get(W + "id")) for e in self.root.iter() if e.get(W + "id") and e.get(W + "id").isdigit()]
        self.next_id = max(ids + [0]) + 1000
        self.log = []

    def nid(self):
        self.next_id += 1
        return str(self.next_id)

    def save(self):
        self.tree.write(self.path, xml_declaration=True, encoding="UTF-8", standalone=True)

    # ---------------- text model ----------------
    def paras(self):
        return list(self.body.iter(W + "p"))

    @staticmethod
    def _deleted(el):
        return any(a.tag in (W + "del", W + "moveFrom") for a in el.iterancestors())

    def segments(self, p):
        """[(w:t element, start, end)] of visible text in paragraph p (nested paragraphs excluded)."""
        segs, pos = [], 0
        for t in p.iter(W + "t"):
            # skip text belonging to a nested paragraph (text boxes)
            if next(t.iterancestors(W + "p")) is not p:
                continue
            if self._deleted(t):
                continue
            s = t.text or ""
            segs.append((t, pos, pos + len(s)))
            pos += len(s)
        return segs

    def text(self, p):
        return "".join(t.text or "" for t, _, _ in self.segments(p))

    def find_para(self, must_contain, nth=None):
        hits = [p for p in self.paras() if must_contain in self.text(p)]
        if nth is None:
            assert len(hits) == 1, (must_contain[:60], len(hits))
            return hits[0]
        return hits[nth]

    # ---------------- structural helpers ----------------
    def _split_run(self, t, off):
        """Split the run holding w:t `t` so that a run boundary falls at character `off` of t.text.
        Returns (left_t, right_t) w:t elements (either may be None when off is 0 / len)."""
        s = t.text or ""
        if off <= 0:
            return None, t
        if off >= len(s):
            return t, None
        r = t.getparent()
        assert r.tag == W + "r", r.tag
        children = list(r)
        k = children.index(t)
        r2 = copy.deepcopy(r)
        # left run keeps children [..k] with truncated text; right run keeps [k..] with the remainder
        for c in list(r)[k + 1:]:
            r.remove(c)
        t.text = s[:off]
        t.set(XMLSPACE, "preserve")
        c2 = list(r2)
        for c in c2[:k]:
            if c.tag != W + "rPr":
                r2.remove(c)
        t2 = [c for c in r2 if c.tag == W + "t"][0]
        t2.text = s[off:]
        t2.set(XMLSPACE, "preserve")
        r.addnext(r2)
        return t, t2

    def _isolate_in_ins(self, runs):
        """runs: consecutive sibling runs sharing a w:ins parent.  Split that w:ins so they sit alone in one part."""
        ins = runs[0].getparent()
        if ins.tag != W + "ins":
            return None
        kids = list(ins)
        i0, i1 = kids.index(runs[0]), kids.index(runs[-1])
        before, mid, after = kids[:i0], kids[i0:i1 + 1], kids[i1 + 1:]
        if before:
            b = etree.Element(W + "ins", ins.attrib)
            b.set(W + "id", self.nid())
            for c in before:
                b.append(c)
            ins.addprevious(b)
        if after:
            a = etree.Element(W + "ins", ins.attrib)
            a.set(W + "id", self.nid())
            for c in after:
                a.append(c)
            ins.addnext(a)
        return ins

    def _make_run(self, rpr, text, italic=None):
        r = etree.Element(W + "r")
        if rpr is not None:
            rp = copy.deepcopy(rpr)
            for bad in rp.findall(W + "ins") + rp.findall(W + "del") + rp.findall(W + "rPrChange"):
                rp.remove(bad)
            if italic is not None:
                for tag in ("i", "iCs"):
                    for e in rp.findall(W + tag):
                        rp.remove(e)
                if italic:
                    # w:i must follow rStyle/rFonts/b/bCs in CT_RPr order
                    anchor = None
                    for tag in ("rStyle", "rFonts", "b", "bCs"):
                        e = rp.find(W + tag)
                        if e is not None:
                            anchor = e
                    ie = etree.Element(W + "i")
                    ice = etree.Element(W + "iCs")
                    if anchor is None:
                        rp.insert(0, ie)
                    else:
                        anchor.addnext(ie)
                    ie.addnext(ice)
            if len(rp):
                r.append(rp)
        elif italic:
            rp = etree.SubElement(r, W + "rPr")
            etree.SubElement(rp, W + "i")
            etree.SubElement(rp, W + "iCs")
        t = etree.SubElement(r, W + "t")
        t.text = text
        t.set(XMLSPACE, "preserve")
        return r

    def _new_runs(self, rpr, new, italic=None):
        """new: plain string, or with *...* marking italic spans (italics toggled relative to rpr).
        italic=True/False forces the italic state of a plain string."""
        if "*" not in new:
            return [self._make_run(rpr, new, italic=italic)] if new else []
        out = []
        parts = re.split(r"(\*[^*]+\*)", new)
        for part in parts:
            if not part:
                continue
            if part.startswith("*") and part.endswith("*"):
                out.append(self._make_run(rpr, part[1:-1], italic=True))
            else:
                out.append(self._make_run(rpr, part, italic=False))
        return out

    # ---------------- tracked replace ----------------
    def replace(self, p, old, new, comment=None, occurrence=None, label="", italic=None):
        txt = self.text(p)
        n = txt.count(old)
        if occurrence is None:
            assert n == 1, ("match count", old[:60], n)
            s = txt.index(old)
        else:
            s = [m.start() for m in re.finditer(re.escape(old), txt)][occurrence]
        e = s + len(old)
        # split runs at s and e
        segs = self.segments(p)
        for pos in (e, s):  # split the end first so start offsets stay valid
            for t, a, b in self.segments(p):
                if a < pos < b:
                    self._split_run(t, pos - a)
                    break
        segs = [(t, a, b) for t, a, b in self.segments(p) if a >= s and b <= e and b > a]
        assert "".join(t.text for t, _, _ in segs) == old, ("reassembly", old[:60])
        runs = []
        for t, _, _ in segs:
            r = t.getparent()
            assert r.tag == W + "r"
            if r not in runs:
                runs.append(r)
        # a matched run must hold only the matched text (no stray tabs etc.)
        for r in runs:
            for c in r:
                if c.tag not in (W + "rPr", W + "t", W + "lastRenderedPageBreak"):
                    raise ValueError(f"matched run holds {c.tag}; handle manually: {old[:60]}")
        rpr = runs[0].find(W + "rPr")
        # group consecutive siblings with the same parent
        groups = []
        for r in runs:
            if groups and groups[-1][-1].getnext() is r:
                groups[-1].append(r)
            else:
                groups.append([r])
        dels = []
        for g in groups:
            par = g[0].getparent()
            assert par.tag in (W + "p", W + "ins"), par.tag
            if par.tag == W + "ins":
                self._isolate_in_ins(g)
            d = etree.Element(W + "del")
            d.set(W + "id", self.nid()); d.set(W + "author", AUTHOR); d.set(W + "date", self.date)
            g[0].addprevious(d)
            for r in g:
                for t in r.findall(W + "t"):
                    t.tag = W + "delText"
                    t.set(XMLSPACE, "preserve")
                d.append(r)
            dels.append(d)
        # anchor elements at paragraph level
        first = dels[0] if dels[0].getparent() is p else dels[0].getparent()
        last = dels[-1] if dels[-1].getparent() is p else dels[-1].getparent()
        assert first.getparent() is p and last.getparent() is p, "nested container not handled"
        ins = None
        newruns = self._new_runs(rpr, new, italic=italic)
        if newruns:
            ins = etree.Element(W + "ins")
            ins.set(W + "id", self.nid()); ins.set(W + "author", AUTHOR); ins.set(W + "date", self.date)
            for r in newruns:
                ins.append(r)
            last.addnext(ins)
        end_anchor = ins if ins is not None else last
        cid = None
        if comment:
            cid = self.comment(first, end_anchor, comment)
        self.log.append(dict(label=label, old=old, new=new, comment=cid))
        return ins

    def insert_after(self, p, anchor_text, new, comment=None, label=""):
        """Tracked insertion of `new` right after the visible text `anchor_text` (unique) in paragraph p."""
        txt = self.text(p)
        assert txt.count(anchor_text) == 1, ("anchor", anchor_text[:60], txt.count(anchor_text))
        pos = txt.index(anchor_text) + len(anchor_text)
        target = None
        for t, a, b in self.segments(p):
            if a < pos < b:
                self._split_run(t, pos - a)
                break
        for t, a, b in self.segments(p):
            if b == pos:
                target = t
        assert target is not None
        r = target.getparent()
        rpr = r.find(W + "rPr")
        par = r.getparent()
        if par.tag == W + "ins":
            # split so that r ends its w:ins part; put the new w:ins after that part
            kids = list(par)
            after = kids[kids.index(r) + 1:]
            if after:
                a = etree.Element(W + "ins", par.attrib); a.set(W + "id", self.nid())
                for c in after:
                    a.append(c)
                par.addnext(a)
            anchor = par
        else:
            assert par is p, par.tag
            anchor = r
        ins = etree.Element(W + "ins")
        ins.set(W + "id", self.nid()); ins.set(W + "author", AUTHOR); ins.set(W + "date", self.date)
        for nr in self._new_runs(rpr, new):
            ins.append(nr)
        anchor.addnext(ins)
        cid = self.comment(ins, ins, comment) if comment else None
        self.log.append(dict(label=label, old="", new=new, comment=cid))
        return ins

    def comment(self, first, last, text):
        """Anchor a new comment from before `first` to after `last` (both direct children of the same w:p)."""
        cid, _, _ = add_comment(self.dir, text, author=AUTHOR, initials=INITIALS)
        cs = etree.Element(W + "commentRangeStart"); cs.set(W + "id", str(cid))
        ce = etree.Element(W + "commentRangeEnd"); ce.set(W + "id", str(cid))
        first.addprevious(cs)
        last.addnext(ce)
        rr = etree.Element(W + "r")
        rp = etree.SubElement(rr, W + "rPr")
        rs = etree.SubElement(rp, W + "rStyle"); rs.set(W + "val", "CommentReference")
        cr = etree.SubElement(rr, W + "commentReference"); cr.set(W + "id", str(cid))
        ce.addnext(rr)
        return cid

    def comment_on_para(self, p, text):
        """Comment spanning the whole paragraph (no text change)."""
        kids = [c for c in p if c.tag != W + "pPr"]
        return self.comment(kids[0], kids[-1], text)

    # ---------------- images ----------------
    def replace_image(self, old_media_name, new_png, comment=None, label=""):
        rels_path = os.path.join(self.dir, "word", "_rels", "document.xml.rels")
        rels = etree.parse(rels_path)
        rid_old = [r.get("Id") for r in rels.getroot() if r.get("Target") == f"media/{old_media_name}"]
        assert len(rid_old) == 1, (old_media_name, rid_old)
        rid_old = rid_old[0]
        blips = [b for b in self.root.iter("{http://schemas.openxmlformats.org/drawingml/2006/main}blip")
                 if b.get("{%s}embed" % RNS) == rid_old]
        assert len(blips) == 1, ("blip uses", old_media_name, len(blips))
        run = next(a for a in blips[0].iterancestors(W + "r"))
        p = next(run.iterancestors(W + "p"))
        # new media part
        media_dir = os.path.join(self.dir, "word", "media")
        k = 1
        while os.path.exists(os.path.join(media_dir, f"image_r232_{k}.png")):
            k += 1
        new_name = f"image_r232_{k}.png"
        data = open(new_png, "rb").read()
        open(os.path.join(media_dir, new_name), "wb").write(data)
        ids = [int(r.get("Id")[3:]) for r in rels.getroot() if r.get("Id", "").startswith("rId") and r.get("Id")[3:].isdigit()]
        rid_new = f"rId{max(ids) + 1}"
        rel = etree.SubElement(rels.getroot(), "{http://schemas.openxmlformats.org/package/2006/relationships}Relationship")
        rel.set("Id", rid_new)
        rel.set("Type", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image")
        rel.set("Target", f"media/{new_name}")
        rels.write(rels_path, xml_declaration=True, encoding="UTF-8", standalone=True)
        # png size -> aspect
        w_px, h_px = struct.unpack(">II", data[16:24])
        new_run = copy.deepcopy(run)
        A = "{http://schemas.openxmlformats.org/drawingml/2006/main}"
        WP = "{http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing}"
        for b in new_run.iter(A + "blip"):
            b.set("{%s}embed" % RNS, rid_new)
        docpr_ids = [int(e.get("id")) for e in self.root.iter(WP + "docPr") if (e.get("id") or "").isdigit()]
        for dp in new_run.iter(WP + "docPr"):
            dp.set("id", str(max(docpr_ids) + 1))
            dp.set("name", f"Picture R232 {k}")
        for ext in list(new_run.iter(WP + "extent")) + list(new_run.iter(A + "ext")):
            if ext.get("cx") and ext.get("cy"):
                cx = int(ext.get("cx"))
                ext.set("cy", str(round(cx * h_px / w_px)))
        for pr in new_run.iter("{http://schemas.openxmlformats.org/drawingml/2006/picture}cNvPr"):
            pr.set("name", new_name)
        # delete old run, insert new
        par = run.getparent()
        if par.tag == W + "ins":
            self._isolate_in_ins([run])
        d = etree.Element(W + "del")
        d.set(W + "id", self.nid()); d.set(W + "author", AUTHOR); d.set(W + "date", self.date)
        run.addprevious(d); d.append(run)
        outer = d if d.getparent() is p else d.getparent()
        assert outer.getparent() is p
        ins = etree.Element(W + "ins")
        ins.set(W + "id", self.nid()); ins.set(W + "author", AUTHOR); ins.set(W + "date", self.date)
        ins.append(new_run)
        outer.addnext(ins)
        cid = self.comment(outer, ins, comment) if comment else None
        self.log.append(dict(label=label, old=f"image {old_media_name}", new=f"image {new_name} ({w_px}x{h_px}) from {new_png}", comment=cid))
        return p


def fix_settings_track_revisions(unpacked):
    """Replace the invalid <w:trackChanges/> (written by an earlier helper at the top of settings.xml) by a
    schema-valid <w:trackRevisions/> in CT_Settings order.  Content is unaffected."""
    path = os.path.join(unpacked, "word", "settings.xml")
    tree = etree.parse(path); root = tree.getroot()
    bad = root.findall(W + "trackChanges")
    for b in bad:
        root.remove(b)
    if root.find(W + "trackRevisions") is None:
        before = ["writeProtection", "view", "zoom", "removePersonalInformation", "removeDateAndTime",
                  "doNotDisplayPageBoundaries", "displayBackgroundShape", "printPostScriptOverText",
                  "printFractionalCharacterWidth", "printFormsData", "embedTrueTypeFonts", "embedSystemFonts",
                  "saveSubsetFonts", "saveFormsData", "mirrorMargins", "alignBordersAndEdges",
                  "bordersDoNotSurroundHeader", "bordersDoNotSurroundFooter", "gutterAtTop", "hideSpellingErrors",
                  "hideGrammaticalErrors", "activeWritingStyle", "proofState", "formsDesign", "attachedTemplate",
                  "linkStyles", "stylePaneFormatFilter", "stylePaneSortMethod", "documentType", "mailMerge",
                  "revisionView"]
        last = None
        for c in root:
            if etree.QName(c).localname in before:
                last = c
        tr = etree.Element(W + "trackRevisions")
        if last is None:
            root.insert(0, tr)
        else:
            last.addnext(tr)
    tree.write(path, xml_declaration=True, encoding="UTF-8", standalone=True)
    return len(bad)
