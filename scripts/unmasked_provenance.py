#!/usr/bin/env python3
"""
Helpers for predict_unmasked (rules/augustus_predict/predict_unmasked.smk).

prefix      Prefix gene and transcript IDs of an AUGUSTUS GTF, so that the
            unmasked predictions cannot collide with the masked ones
            (both use g1, g1.t1, ...) when TSEBRA reads them together.

              unmasked_provenance.py prefix --gtf IN --prefix unm_ --out OUT

provenance  Label each transcript of the final gene set by the source gene
            sets that contain a transcript with the identical coding
            structure (sequence, strand, CDS plus stop codon intervals). Independent of
            the IDs, which TSEBRA and rename_gtf.py replace.

              unmasked_provenance.py provenance --final braker.gtf \\
                  --source masked_augustus=augustus.hints.fixed.gtf \\
                  --source unmasked_augustus=unmasked/augustus.hints.fixed.prefixed.gtf \\
                  --source genemark=GeneMark-ETP/genemark.gtf \\
                  --source training=GeneMark-ETP/training.gtf \\
                  --out provenance.tsv --summary provenance_summary.txt

            class (last column):
              masked         found in a gene set predicted with repeats
                             masked (masked AUGUSTUS, GeneMark, training)
              unmasked_only  found only in the unmasked AUGUSTUS set, i.e.
                             the model exists because of the unmasked run
              unmatched      no source has the identical CDS structure
                             (changed after TSEBRA, e.g. by CDS
                             normalisation)
"""
import argparse
import re
import sys
from collections import Counter, defaultdict

TX_RE = re.compile(r'transcript_id "([^"]+)"')
GENE_RE = re.compile(r'gene_id "([^"]+)"')
MASKED_SOURCES = ("masked_augustus", "genemark", "training")


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
    """transcript_id -> (seqid, strand, coding intervals), and transcript -> gene.

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


def cmd_provenance(a):
    final, gene_of = cds_structures(a.final)
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
            hits = [n for n, ks in source_keys.items() if key in ks]
            if any(h in MASKED_SOURCES for h in hits):
                cls = "masked"
            elif "unmasked_augustus" in hits:
                cls = "unmasked_only"
            else:
                cls = "unmatched"
            counts[cls] += 1
            out.write(f"{gene_of.get(tx, '')}\t{tx}\t{key[0]}\t{key[1]}\t{key[2][0][0]}\t"
                      f"{key[2][-1][1]}\t{len(key[2])}\t{','.join(sorted(hits)) or '-'}\t{cls}\n")

    lines = [f"transcripts\t{len(final)}"] + [f"{k}\t{counts[k]}" for k in ("masked", "unmasked_only", "unmatched")]
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
    p = sub.add_parser("provenance")
    p.add_argument("--final", required=True)
    p.add_argument("--source", action="append", required=True, help="NAME=GTF (repeatable)")
    p.add_argument("--out", required=True)
    p.add_argument("--summary", required=True)
    a = ap.parse_args()
    cmd_prefix(a) if a.cmd == "prefix" else cmd_provenance(a)


if __name__ == "__main__":
    main()
