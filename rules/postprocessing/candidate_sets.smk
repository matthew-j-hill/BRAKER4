"""
Extra TSEBRA candidate gene sets and the provenance of the final transcripts.

Candidate sets come from predict_unmasked = 1 (AUGUSTUS on the unmasked
genome, rules/augustus_predict/predict_unmasked.smk) and from
[paths] extra_candidate_gtf = name:path[,name:path...], gene models from any
other predictor (GTF or GFF3). Each extra set is converted to an AUGUSTUS-style
GTF with IDs prefixed '<name>_' (so they cannot collide with AUGUSTUS or
GeneMark IDs) and given to TSEBRA as a --gtf candidate set, filtered against
the hints like AUGUSTUS and GeneMark (run_tsebra.smk, candidate_gtfs).

provenance/provenance.tsv labels each transcript of the final braker.gtf by
the gene sets with the identical coding structure (CDS plus stop codon):
'masked' if a set predicted on the soft-masked genome has it (masked AUGUSTUS,
GeneMark, training genes), otherwise the extra set(s) that have it (e.g.
'unmasked_augustus', 'vipsania'), or 'unmatched'. collect_results copies it to
results/provenance/.
"""


if config['extra_candidates']:

    rule prepare_candidate_set:
        """
        Convert an extra candidate gene set (GTF or GFF3) to an AUGUSTUS-style
        GTF for TSEBRA: gene/transcript lines, CDS, introns, start/stop codons,
        IDs prefixed '<name>_'.
        """
        input:
            lambda w: config['extra_candidates'][w.cand]
        output:
            gtf = "output/{sample}/candidates/{cand}.gtf"
        wildcard_constraints:
            cand = "|".join(config['extra_candidates'])
        threads: 1
        resources:
            mem_mb=mem_single_thread,
            runtime=int(config['slurm_args']['max_runtime'])
        container:
            BRAKER3_CONTAINER
        shell:
            r"""
            export PATH=/opt/conda/bin:$PATH
            export PYTHONNOUSERSITE=1
            set -euo pipefail
            python3 {script_dir}/unmasked_provenance.py candidate \
                --in {input} --prefix {wildcards.cand}_ --source {wildcards.cand} --out {output.gtf}
            """


def _provenance_sources(wildcards):
    """Gene sets the final transcripts are matched against. GeneMark and the
    training genes come from the same mode-aware inputs as TSEBRA's."""
    sample = wildcards.sample
    t = get_tsebra_inputs(wildcards)
    sources = {
        "final": f"output/{sample}/braker.gtf",
        "masked_augustus": f"output/{sample}/augustus.hints.fixed.gtf",
        "genemark": t["genemark_gtf"],
        "training": t["training_gtf"],
    }
    for label, gtf in candidate_gtfs(sample):
        sources["unmasked_augustus" if label == "unmasked" else label] = gtf
    # decisions of the automatic mode fallback, as for run_tsebra
    sources["mode_decisions"] = mode_decision_files(wildcards)
    return sources


def _provenance_source_args(wildcards, input):
    return " ".join(f"--source {k}={input[k]}" for k in input.keys()
                    if k not in ("final", "mode_decisions"))


rule candidate_provenance:
    """
    Label each transcript of the final gene set by the gene sets with the
    identical coding structure.
    """
    input:
        unpack(_provenance_sources)
    output:
        tsv = "output/{sample}/provenance/provenance.tsv",
        summary = "output/{sample}/provenance/provenance_summary.txt"
    params:
        sources = _provenance_source_args
    threads: 1
    resources:
        mem_mb=mem_single_thread,
        runtime=int(config['slurm_args']['max_runtime'])
    container:
        BRAKER3_CONTAINER
    shell:
        r"""
        export PATH=/opt/conda/bin:$PATH
        export PYTHONNOUSERSITE=1
        set -euo pipefail
        python3 {script_dir}/unmasked_provenance.py provenance \
            --final {input.final} {params.sources} \
            --masked masked_augustus,genemark,training \
            --out {output.tsv} --summary {output.summary}
        """
