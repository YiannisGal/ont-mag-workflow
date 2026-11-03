# ---------------------------------------------------------------------
# GTDB-Tk classification (incl. RED) and skani ANI
# ---------------------------------------------------------------------
REFS = config.get("reference_genomes", "")


rule gtdbtk_classify:
    input:
        f"{OUT}/mags/fasta",
    output:
        f"{OUT}/taxonomy/gtdbtk_summary.tsv",
    params:
        outdir=f"{OUT}/taxonomy/gtdbtk",
        db=DB["gtdbtk"],
        ani=(f"--mash_db {DB['gtdbtk_mash_db']}" if is_set(DB.get("gtdbtk_mash_db")) else "--skip_ani_screen"),
    log:
        f"{OUT}/logs/taxonomy/gtdbtk.log",
    threads: T["big"]
    resources:
        mem_mb=MEM["gtdbtk"],
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        """
        export GTDBTK_DATA_PATH={params.db}
        gtdbtk classify_wf --genome_dir {input} -x fa --out_dir {params.outdir} \
            --cpus {threads} {params.ani} > {log} 2>&1
        # merge bacterial + archaeal summaries (either may be absent)
        first=1
        for f in {params.outdir}/gtdbtk.bac120.summary.tsv {params.outdir}/gtdbtk.ar53.summary.tsv; do
            [ -s "$f" ] || continue
            if [ $first -eq 1 ]; then cat "$f"; first=0; else tail -n +2 "$f"; fi
        done > {output}
        """


rule skani_mags:
    """All-vs-all ANI among MAGs. skani omits pairs below its screen (~80% ANI):
    treat missing pairs as 'below the lowest value reported in the run', not zero."""
    input:
        f"{OUT}/mags/fasta",
    output:
        f"{OUT}/taxonomy/skani_mags.tsv",
    threads: T["small"]
    conda:
        "../envs/skani.yaml"
    shell:
        "skani triangle -t {threads} --sparse {input}/*.fa > {output}"


rule skani_refs:
    input:
        f"{OUT}/mags/fasta",
    output:
        f"{OUT}/taxonomy/skani_mags_vs_refs.tsv",
    params:
        refs=REFS,
    threads: T["small"]
    conda:
        "../envs/skani.yaml"
    shell:
        """
        shopt -s nullglob
        skani dist -t {threads} -q {input}/*.fa -r {params.refs}/*.fa {params.refs}/*.fna {params.refs}/*.fasta > {output}
        """
