#!/usr/bin/env python
"""Local DAS Tool 1.1.7 result -> the exact Workspace.save_objects params for one KBaseMetagenomes.BinnedContigs.

Builds (does NOT save) the object that MetagenomeUtils/extract_bins_as_assemblies (beta 54c6bc35) consumes:
binned_contigs_to_file() reads only data.assembly_ref, bins[].bid and the keys of bins[].contigs, writes each bin as
'>contig_id\\nSEQUENCE.upper()\\n' from the Assembly's blob FASTA, and AssemblyUtil.save_assembly_from_fasta names the
new Assembly bid + assembly_suffix. Contig gc/len are taken from the Assembly object (as KBase binners do:
MetagenomeFileUtils._generate_contigs), bin gc / sum_contig_len are computed from the sequences, cov = 0.0 (as
kb_das_tool's own output; kb_Msuite's withFilter path reads bins[].cov).
Also predicts, per bin, the Assembly that extract will save (name, num_contigs, dna_size, AssemblyUtil md5 =
md5(','.join(sorted(md5(seq.upper()) per contig)))) so the driver can verify extract's output without downloading.
The data object is validated against the Workspace's own JSON schema for KBaseMetagenomes.BinnedContigs-1.0.
Usage: python -I dastool_to_binnedcontigs.py --fasta A.fna --assembly-upa 275645/X/1 --dastool-prefix out/S
          --sample S --binner metabat2=REF --binner maxbin2=REF --binner concoct=REF --out S.save_objects.json
"""
import argparse, hashlib, json, os, re, socket, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kbrpc import rpc

WS_ID, TYPE = 275645, "KBaseMetagenomes.BinnedContigs-1.0"
NAME_OK = re.compile(r"^[A-Za-z0-9|._-]+$")          # Workspace object-name alphabet


def read_fasta(path):
    seqs, cur, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf).upper()
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf).upper()
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--assembly-upa", required=True)
    ap.add_argument("--dastool-prefix", required=True, help="DAS Tool -o value (reads <prefix>_DASTool_contig2bin.tsv)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--binner", action="append", default=[], help="label=BinnedContigs UPA (provenance inputs)")
    ap.add_argument("--ws-id", type=int, default=WS_ID)
    ap.add_argument("--assembly-suffix", default=None, help="extract's assembly_suffix (default __<sample>)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    suffix = a.assembly_suffix if a.assembly_suffix is not None else f"__{a.sample}"
    c2b = a.dastool_prefix + "_DASTool_contig2bin.tsv"
    summ = a.dastool_prefix + "_DASTool_summary.tsv"
    if not (os.path.isfile(c2b) and os.path.getsize(c2b) and os.path.isfile(summ)):
        sys.exit(f"no DAS Tool bins at {a.dastool_prefix} (no score >= threshold?) -> mark the stage 'empty'")
    seqs = read_fasta(a.fasta)
    asm = rpc("ws", "Workspace.get_objects2", [{"objects": [{"ref": a.assembly_upa, "included": [
        "/contigs/*/gc_content", "/contigs/*/length", "/md5", "/num_contigs"]}]}])[0]["data"][0]
    if f"{asm['info'][6]}/{asm['info'][0]}/{asm['info'][4]}" != a.assembly_upa:
        sys.exit("assembly UPA must be absolute (wsid/objid/ver)")
    oc = asm["data"]["contigs"]
    bins, order = {}, []
    with open(c2b) as fh:
        for line in fh:
            cid, lab = line.rstrip("\n").split("\t")
            if lab == "unbinned":
                continue
            if lab not in bins:
                bins[lab] = []
                order.append(lab)
            bins[lab].append(cid)
    with open(summ) as fh:                      # summary rows are in DAS Tool selection order (best score first)
        hdr = fh.readline().rstrip("\n").split("\t")
        summary = {r[0]: dict(zip(hdr, r)) for r in (l.rstrip("\n").split("\t") for l in fh if l.strip())}
    assert set(summary) == set(bins), ("summary/contig2bin bins differ", set(summary) ^ set(bins))
    seen, out_bins, predicted = set(), [], []
    for lab in sorted(bins, key=lambda b: -float(summary[b]["bin_score"])):
        bid = f"{lab}.fasta"
        if "," in bid or not NAME_OK.match(bid + suffix):
            sys.exit(f"bin label {lab!r} is not usable as a bid / Workspace object name")
        contigs, gcn, tot, md5s = {}, 0, 0, []
        for cid in bins[lab]:
            if cid in seen:
                sys.exit(f"contig {cid} in two DAS Tool bins")
            seen.add(cid)
            s = seqs[cid]
            if oc[cid]["length"] != len(s):
                sys.exit(f"{cid}: Assembly length {oc[cid]['length']} != FASTA {len(s)}")
            contigs[cid] = {"gc": oc[cid]["gc_content"], "len": oc[cid]["length"]}
            gcn += s.count("G") + s.count("C")
            tot += len(s)
            md5s.append(hashlib.md5(s.encode()).hexdigest())
        out_bins.append({"bid": bid, "contigs": contigs, "n_contigs": len(contigs), "gc": round(gcn / tot, 5),
                         "sum_contig_len": tot, "cov": 0.0})
        predicted.append({"bid": bid, "dastool_label": lab, "assembly_name": bid + suffix, "num_contigs": len(contigs),
                          "dna_size": tot, "md5": hashlib.md5(",".join(sorted(md5s)).encode()).hexdigest(),
                          "SCG_completeness": summary[lab]["SCG_completeness"],
                          "SCG_redundancy": summary[lab]["SCG_redundancy"], "bin_score": summary[lab]["bin_score"]})
    data = {"assembly_ref": a.assembly_upa, "bins": out_bins, "total_contig_len": sum(b["sum_contig_len"] for b in out_bins)}
    ti = rpc("ws", "Workspace.get_type_info", [TYPE])[0]
    import jsonschema
    jsonschema.validate(data, json.loads(ti["json_schema"]))
    script = os.path.abspath(__file__)
    prov = [{"time": time.strftime("%Y-%m-%dT%H:%M:%S+0000", time.gmtime()),
             "script": "pipeline.py dastool (local DAS Tool 1.1.7)",
             "script_ver": hashlib.sha256(open(script, "rb").read()).hexdigest()[:12],
             "script_command_line": " ".join(["python", "-I"] + sys.argv),
             "description": (f"DAS Tool 1.1.7 run locally on {socket.gethostname()} over "
                             + ", ".join(b.split('=')[0] for b in a.binner)
                             + " BinnedContigs of the same Assembly: --search_engine diamond --score_threshold 0.5 "
                               "--duplicate_penalty 0.6 --megabin_penalty 0.5; bids = DAS Tool bin labels + '.fasta'"),
             "input_ws_objects": [a.assembly_upa] + [b.split("=", 1)[1] for b in a.binner]}]
    params = [{"id": a.ws_id, "objects": [{"type": TYPE, "name": f"{a.sample}.DASTool.BinnedContigs", "data": data,
                                           "meta": {"DAS_Tool": "1.1.7 (local)"}, "provenance": prov}]}]
    json.dump({"save_objects_params": params, "predicted_assemblies": predicted,
               "extract_params": [{"workspace_name": "<WS_NAME>", "binned_contig_obj_ref": "<saved BinnedContigs UPA>",
                                   "extracted_assemblies": ",".join(b["bid"] for b in out_bins),
                                   "assembly_suffix": suffix, "assembly_set_name": f"{a.sample}.DASTool_bins.AssemblySet"}]},
              open(a.out, "w"), indent=1)
    print(json.dumps({"bins": len(out_bins), "binned_bp": data["total_contig_len"], "schema": "valid", "out": a.out,
                      "first_bin": {k: v for k, v in out_bins[0].items() if k != "contigs"},
                      "predicted_first": predicted[0]}, indent=1))


if __name__ == "__main__":
    main()
