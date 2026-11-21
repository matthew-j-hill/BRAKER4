"""
AUGUSTUS prediction on the unmasked genome ([PARAMS] predict_unmasked = 1).

Training is unchanged: GeneMark and AUGUSTUS are trained on the soft-masked
genome as usual. AUGUSTUS then predicts a second time with the same
parameters, hints and extrinsic config, but with --softmasking=0, so lowercase
(repeat) sequence is not penalised. Genes truncated or missed because an exon
overlaps masked sequence can then be predicted.

The unmasked gene set enters TSEBRA as one more CANDIDATE set (--gtf),
filtered against the hints like masked AUGUSTUS and GeneMark: an unmasked
model is kept only where the evidence supports it at least as well as the
alternatives. Only ETP-style modes (etp, isoseq, dual) filter AUGUSTUS in
TSEBRA; in ES/EP/ET the unmasked set is not used (run_tsebra.smk,
unmasked_candidate_gtf).

All files are under output/{sample}/unmasked/. The existing rules are reused
with `use rule ... with:`, which replaces the directives listed here; input,
output and params are given in full so the result does not depend on whether
Snakemake merges or replaces them.

Provenance: unmasked/provenance.tsv labels each transcript of the final
braker.gtf as masked / unmasked_only / unmatched by identical CDS structure
(scripts/unmasked_provenance.py); collect_results copies it to results/.
"""


use rule run_augustus_hints as run_augustus_hints_unmasked with:
    input:
        genome = lambda w: get_masked_genome(w.sample),
        hintsfile = "output/{sample}/hintsfile.gff",
        optimize_log = "output/{sample}/optimize_augustus.log",
        mode_decisions = mode_decision_files
    output:
        job_lst = "output/{sample}/unmasked/augustus_hints.job.lst",
        augustus_gff = "output/{sample}/unmasked/augustus.hints.gff",
        augustus_gtf = "output/{sample}/unmasked/augustus.hints.gtf"
    benchmark:
        "benchmarks/{sample}/run_augustus_hints_unmasked/run_augustus_hints_unmasked.txt"
    params:
        species_name = lambda w: get_species_name(w),
        output_dir = lambda w: f"output/{w.sample}/unmasked",
        aug_config = augustus_config_path,
        extrinsic_cfg = lambda w: get_extrinsic_cfg(w),
        chunksize = config['augustus_chunksize'],
        overlap = config['augustus_overlap'],
        use_dev_shm = config['use_dev_shm'],
        allow_hinted_splicesites = config.get('allow_hinted_splicesites', 'gcag,atac'),
        tmp_root = lambda w: "/dev/shm" if config['use_dev_shm'] else TMP_ROOT,
        softmasking = 0


use rule get_anno_fasta as get_anno_fasta_unmasked with:
    input:
        gtf = "output/{sample}/unmasked/augustus.hints.gtf",
        genome = lambda w: os.path.join(get_braker_dir(w), "genome.fa"),
        mode_decisions = mode_decision_files
    output:
        codingseq = "output/{sample}/unmasked/augustus.hints.codingseq",
        aa = "output/{sample}/unmasked/augustus.hints.aa",
        bad_genes = "output/{sample}/unmasked/bad_genes.lst"
    benchmark:
        "benchmarks/{sample}/get_anno_fasta_unmasked/get_anno_fasta_unmasked.txt"
    params:
        output_dir = lambda w: f"output/{w.sample}/unmasked",
        translation_table = config.get("translation_table", 1)


use rule fix_in_frame_stop_codons as fix_in_frame_stop_codons_unmasked with:
    input:
        gtf = "output/{sample}/unmasked/augustus.hints.gtf",
        bad_genes = "output/{sample}/unmasked/bad_genes.lst",
        genome = lambda w: os.path.join(get_braker_dir(w), "genome.fa"),
        hintsfile = get_augustus_hintsfile,
        species_params = "augustus_config/species/{sample}_rescuedB3/{sample}_rescuedB3_parameters.cfg.backup",
        mode_decisions = mode_decision_files
    output:
        gtf_fixed = "output/{sample}/unmasked/augustus.hints.fixed.gtf",
        codingseq_fixed = "output/{sample}/unmasked/augustus.hints.fixed.codingseq",
        aa_fixed = "output/{sample}/unmasked/augustus.hints.fixed.aa",
        fix_log = "output/{sample}/unmasked/fix_in_frame_stop_codons.log"
    benchmark:
        "benchmarks/{sample}/fix_in_frame_stop_codons_unmasked/fix_in_frame_stop_codons_unmasked.txt"
    params:
        species_name = lambda w: get_species_name(w),
        output_dir = lambda w: f"output/{w.sample}/unmasked",
        translation_table = config.get("translation_table", 1),
        aug_config = augustus_config_path,
        extrinsic_cfg = lambda w: get_extrinsic_cfg(w),
        fix_stem = lambda w: f"output/{w.sample}/unmasked/augustus.hints_fix_ifs_",
        allow_hinted_splicesites = config.get('allow_hinted_splicesites', 'gcag,atac'),
        softmask = "off"


rule prefix_unmasked_ids:
    """
    Prefix the unmasked gene and transcript IDs with 'unm_' so they cannot
    collide with the masked AUGUSTUS IDs (both start at g1) inside TSEBRA.
    """
    input:
        gtf = "output/{sample}/unmasked/augustus.hints.fixed.gtf"
    output:
        gtf = "output/{sample}/unmasked/augustus.hints.fixed.prefixed.gtf"
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
        python3 {script_dir}/unmasked_provenance.py prefix \
            --gtf {input.gtf} --prefix unm_ --out {output.gtf}
        """


def _provenance_sources(wildcards):
    """Gene sets the final transcripts are matched against. GeneMark and the
    training genes come from the same mode-aware inputs as TSEBRA's."""
    t = get_tsebra_inputs(wildcards)
    return {
        "final": f"output/{wildcards.sample}/braker.gtf",
        "masked": f"output/{wildcards.sample}/augustus.hints.fixed.gtf",
        "unmasked": f"output/{wildcards.sample}/unmasked/augustus.hints.fixed.prefixed.gtf",
        "genemark": t["genemark_gtf"],
        "training": t["training_gtf"],
        # decisions of the automatic mode fallback, as for run_tsebra
        "mode_decisions": mode_decision_files(wildcards),
    }


rule unmasked_provenance:
    """
    Label each transcript of the final gene set by the gene sets with the
    identical CDS structure: masked (masked AUGUSTUS, GeneMark, training
    genes), unmasked_only, or unmatched.
    """
    input:
        unpack(_provenance_sources)
    output:
        tsv = "output/{sample}/unmasked/provenance.tsv",
        summary = "output/{sample}/unmasked/provenance_summary.txt"
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
            --final {input.final} \
            --source masked_augustus={input.masked} \
            --source unmasked_augustus={input.unmasked} \
            --source genemark={input.genemark} \
            --source training={input.training} \
            --out {output.tsv} --summary {output.summary}
        """
