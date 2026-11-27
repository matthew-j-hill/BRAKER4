#!/usr/bin/env python3
"""
Helpers for extra TSEBRA candidate gene sets (predict_unmasked,
extra_candidate_gtf; rules/augustus_predict/predict_unmasked.smk,
rules/postprocessing/candidate_sets.smk).

prefix      Prefix gene and transcript IDs of an AUGUSTUS GTF, so that the
            unmasked predictions cannot collide with the masked ones
            (both use g1, g1.t1, ...) when TSEBRA reads them together.

              unmasked_provenance.py prefix --gtf IN --prefix unm_ --out OUT

candidate   Convert another predictor's gene models (GTF or GFF3) to an
            AUGUSTUS-style GTF for TSEBRA: gene and transcript lines, CDS,
            introns between CDS, start and stop codons, IDs prefixed. Start
            and stop codons are taken from the input if present, otherwise
            from the CDS ends: by default the CDS is taken to include the stop
            codon (GFF3 convention); --cds-excludes-stop places it after the CDS.

              unmasked_provenance.py candidate --in vipsania.gff3 \\
                  --prefix vipsania_ --source vipsania --out vipsania.gtf

provenance  Label each transcript of the final gene set by the source gene
            sets that contain a transcript with the identical coding
            structure (sequence, strand, CDS plus stop codon intervals).
            Independent of the IDs, which TSEBRA and rename_gtf.py replace.

              unmasked_provenance.py provenance --final braker.gtf \\
                  --source masked_augustus=augustus.hints.fixed.gtf \\
                  --source genemark=GeneMark-ETP/genemark.gtf \\
                  --source training=GeneMark-ETP/training.gtf \\
                  --source unmasked_augustus=unmasked/augustus.hints.fixed.prefixed.gtf \\
                  --source vipsania=candidates/vipsania.gtf \\
                  --masked masked_augustus,genemark,training \\
                  --out provenance.tsv --summary provenance_summary.txt

            class (last column):
              masked         found in a gene set predicted on the soft-masked
                             genome (the --masked sources)
              <set>[+<set>]  found only in the other source(s), e.g.
                             unmasked_augustus, vipsania: the model is in the
                             final set because of that extra set
              unmatched      no source has the identical coding structure
"""
import argparse
import re
import sys
from collections import Counter, defaultdict

TX_RE = re.compile(r'transcript_id "([^"]+)"')
GENE_RE = re.compile(r'gene_id "([^"]+)"')
ATTR_RE = re.compile(r'(?:^|;)\s*([^=;]+)=([^;]*)')


def _merge(intervals):
    """Merge overlapping or adjacent (end + 1 == start) intervals."""
    out = []
    for start, end in sorted(intervals):
        if out and start <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return tuple((a, b) for a, b in out)


def cds_structures(path):
    """GTF: transcript_id -> (seqid, strand, coding intervals), and transcript -> gene.

    Coding intervals are CDS plus stop_codon features, merged. Gene sets differ
    in whether the stop codon is part of the CDS (GeneMark, AUGUSTUS) or only a
    separate stop_codon feature (BRAKER's final braker.gtf after
    normalize_cds trims it); the merged span is the same either way.
    """
    iv = defaultdict(list)
    where = {}
    gene_of = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] not in ("CDS", "stop_codon"):
                continue
            m = TX_RE.search(c[8])
            if not m:
                continue
            tx = m.group(1)
            iv[tx].append((int(c[3]), int(c[4])))
            if c[2] == "CDS":
                where[tx] = (c[0], c[6])
                g = GENE_RE.search(c[8])
                if g:
                    gene_of[tx] = g.group(1)
    keys = {tx: (where[tx][0], where[tx][1], _merge(iv[tx])) for tx in where}
    return keys, gene_of


def cmd_prefix(a):
    n = 0
    with open(a.gtf) as f, open(a.out, "w") as out:
        for line in f:
            if line.startswith("#") or not line.strip():
                out.write(line)
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) >= 9:
                if '"' not in c[8] and c[2] in ("gene", "transcript"):
                    # AUGUSTUS gene/transcript lines carry the bare ID
                    c[8] = a.prefix + c[8].strip()
                else:
                    c[8] = re.sub(r'(gene_id|transcript_id) "', r'\1 "' + a.prefix, c[8])
                n += 1
            out.write("\t".join(c) + "\n")
    print(f"[prefix] {n} feature lines prefixed with '{a.prefix}'", file=sys.stderr)


def read_models(path):
    """Gene models from GTF or GFF3: transcript -> dict(seq, strand, gene,
    cds=[(s, e, phase)], start=[(s, e)], stop=[(s, e)])."""
    models = {}
    tx_gene = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            ftype, attrs = c[2], c[8]
            m = TX_RE.search(attrs)
            if m:                                    # GTF
                txs = [m.group(1)]
                g = GENE_RE.search(attrs)
                if g:
                    tx_gene[txs[0]] = g.group(1)
            else:                                    # GFF3
                kv = dict((k.strip(), v) for k, v in ATTR_RE.findall(attrs))
                if ftype in ("mRNA", "transcript") and "ID" in kv:
                    if "Parent" in kv:
                        tx_gene[kv["ID"]] = kv["Parent"].split(",")[0]
                    continue
                if "Parent" not in kv:
                    continue
                txs = kv["Parent"].split(",")
            if ftype not in ("CDS", "start_codon", "stop_codon"):
                continue
            for tx in txs:
                d = models.setdefault(tx, {"seq": c[0], "strand": c[6], "cds": [], "start": [], "stop": []})
                s, e = int(c[3]), int(c[4])
                if ftype == "CDS":
                    d["cds"].append((s, e, c[7]))
                elif ftype == "start_codon":
                    d["start"].append((s, e))
                else:
                    d["stop"].append((s, e))
    for tx, d in models.items():
        d["gene"] = tx_gene.get(tx, tx)
    return {tx: d for tx, d in models.items() if d["cds"]}


def cmd_candidate(a):
    models = read_models(a.inp)
    genes = defaultdict(list)
    for tx, d in models.items():
        d["cds"].sort()
        genes[d["gene"]].append(tx)
    derived = 0
    lines = []
    order = sorted(genes, key=lambda g: (models[genes[g][0]]["seq"], min(models[t]["cds"][0][0] for t in genes[g])))
    for gene in order:
        txs = genes[gene]
        g = a.prefix + gene
        seq = models[txs[0]]["seq"]
        strand = models[txs[0]]["strand"]
        gstart = min(models[t]["cds"][0][0] for t in txs)
        gend = max(models[t]["cds"][-1][1] for t in txs)
        lines.append((seq, gstart, 0, f"{seq}\t{a.source}\tgene\t{gstart}\t{gend}\t.\t{strand}\t.\t{g}"))
        for tx in txs:
            d = models[tx]
            t = a.prefix + tx
            cds = d["cds"]
            s0, e0 = cds[0][0], cds[-1][1]
            attr = f'transcript_id "{t}"; gene_id "{g}";'
            if a.cds_excludes_stop:
                start = d["start"] or [(s0, s0 + 2) if strand == "+" else (e0 - 2, e0)]
                stop = d["stop"] or [(e0 + 1, e0 + 3) if strand == "+" else (s0 - 3, s0 - 1)]
            else:
                start = d["start"] or [(s0, s0 + 2) if strand == "+" else (e0 - 2, e0)]
                stop = d["stop"] or [(e0 - 2, e0) if strand == "+" else (s0, s0 + 2)]
            if not d["stop"]:
                derived += 1
            ts, te = min(s0, *(x[0] for x in stop)), max(e0, *(x[1] for x in stop))
            lines.append((seq, ts, 1, f"{seq}\t{a.source}\ttranscript\t{ts}\t{te}\t.\t{strand}\t.\t{t}"))
            feats = [("CDS", s, e, ph) for s, e, ph in cds]
            feats += [("intron", x[1] + 1, y[0] - 1, ".") for x, y in zip(cds, cds[1:]) if y[0] - x[1] > 1]
            feats += [("start_codon", s, e, "0") for s, e in start]
            feats += [("stop_codon", s, e, "0") for s, e in stop]
            for ftype, s, e, ph in sorted(feats, key=lambda x: (x[1], x[2])):
                lines.append((seq, s, 2, f"{seq}\t{a.source}\t{ftype}\t{s}\t{e}\t.\t{strand}\t{ph}\t{attr}"))
    # keep each gene block together, genes ordered by position
    with open(a.out, "w") as out:
        for line in (x[3] for x in lines):
            out.write(line + "\n")
    print(f"[candidate] {len(genes)} genes, {len(models)} transcripts from {a.inp}; "
          f"start/stop codons derived from CDS ends for {derived} transcripts", file=sys.stderr)


def cmd_provenance(a):
    final, gene_of = cds_structures(a.final)
    masked = set(a.masked.split(",")) if a.masked else set()
    source_keys = {}
    for spec in a.source:
        name, path = spec.split("=", 1)
        keys, _ = cds_structures(path)
        source_keys.setdefault(name, set()).update(keys.values())

    counts = Counter()
    with open(a.out, "w") as out:
        out.write("gene_id\ttranscript_id\tseqid\tstrand\tcds_start\tcds_end\tn_cds\tsources\tclass\n")
        for tx in sorted(final, key=lambda t: (final[t][0], final[t][2][0][0])):
            key = final[tx]
            hits = sorted(n for n, ks in source_keys.items() if key in ks)
            if any(h in masked for h in hits):
                cls = "masked"
            elif hits:
                cls = "+".join(hits)
            else:
                cls = "unmatched"
            counts[cls] += 1
            out.write(f"{gene_of.get(tx, '')}\t{tx}\t{key[0]}\t{key[1]}\t{key[2][0][0]}\t"
                      f"{key[2][-1][1]}\t{len(key[2])}\t{','.join(hits) or '-'}\t{cls}\n")

    lines = [f"transcripts\t{len(final)}"] + [f"{k}\t{v}" for k, v in counts.most_common()]
    with open(a.summary, "w") as s:
        s.write("\n".join(lines) + "\n")
    print("[provenance] " + "; ".join(lines), file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prefix")
    p.add_argument("--gtf", required=True)
    p.add_argument("--prefix", default="unm_")
    p.add_argument("--out", required=True)
    p = sub.add_parser("candidate")
    p.add_argument("--in", dest="inp", required=True, help="GTF or GFF3")
    p.add_argument("--prefix", required=True)
    p.add_argument("--source", required=True, help="source column (column 2)")
    p.add_argument("--cds-excludes-stop", action="store_true",
                   help="input CDS ends before the stop codon (default: CDS includes it)")
    p.add_argument("--out", required=True)
    p = sub.add_parser("provenance")
    p.add_argument("--final", required=True)
    p.add_argument("--source", action="append", required=True, help="NAME=GTF (repeatable)")
    p.add_argument("--masked", default="masked_augustus,genemark,training",
                   help="comma-separated sources predicted on the soft-masked genome")
    p.add_argument("--out", required=True)
    p.add_argument("--summary", required=True)
    a = ap.parse_args()
    {"prefix": cmd_prefix, "candidate": cmd_candidate, "provenance": cmd_provenance}[a.cmd](a)


if __name__ == "__main__":
    main()
