# ---------------------------------------------------------------------
# Pool DAS_Tool bins -> CheckM2, GUNC, barrnap, tRNAscan-SE ->
# MIMAG classification, filtering, renaming (checkpoint)
# ---------------------------------------------------------------------


rule pool_bins:
    input:
        expand(f"{OUT}/binning/{{asm}}/dastool/{{asm}}_DASTool_bins", asm=ASM),
    output:
        directory(f"{OUT}/bins/candidates"),
    params:
        pairs=[f"{a}:{OUT}/binning/{a}/dastool/{a}_DASTool_bins" for a in ASM],
    shell:
        """
        mkdir -p {output}
        shopt -s nullglob
        for p in {params.pairs}; do
            s=${{p%%:*}}; d=${{p#*:}}
            for f in "$d"/*.fa; do
                b=$(basename "$f" .fa | tr '.' '_')
                cp "$f" {output}/"$s"__"$b".fa
            done
        done
        """


rule checkm2:
    input:
        f"{OUT}/bins/candidates",
    output:
        f"{OUT}/bins/checkm2/quality_report.tsv",
    params:
        outdir=f"{OUT}/bins/checkm2",
        db=DB["checkm2"],
    log:
        f"{OUT}/logs/bins/checkm2.log",
    threads: T["big"]
    conda:
        "../envs/checkm2.yaml"
    shell:
        "checkm2 predict --input {input} -x fa --output-directory {params.outdir} "
        "--threads {threads} --database_path {params.db} --force > {log} 2>&1"


rule gunc:
    input:
        f"{OUT}/bins/candidates",
    output:
        directory(f"{OUT}/bins/gunc"),
    params:
        db=DB["gunc"],
        tmp=f"{OUT}/bins/gunc_tmp",
    log:
        f"{OUT}/logs/bins/gunc.log",
    threads: T["big"]
    conda:
        "../envs/gunc.yaml"
    shell:
        """
        mkdir -p {output} {params.tmp}
        gunc run --input_dir {input} --file_suffix .fa --db_file {params.db} \
            --out_dir {output} --temp_dir {params.tmp} --threads {threads} > {log} 2>&1
        rm -rf {params.tmp}
        """


rule rrna_trna:
    input:
        f"{OUT}/bins/candidates",
    output:
        directory(f"{OUT}/bins/rrna_trna"),
    log:
        f"{OUT}/logs/bins/rrna_trna.log",
    threads: T["small"]
    conda:
        "../envs/rrna_trna.yaml"
    shell:
        """
        mkdir -p {output}
        shopt -s nullglob
        for f in {input}/*.fa; do
            b=$(basename "$f" .fa)
            barrnap --kingdom bac --threads {threads} --quiet "$f" > {output}/"$b".bac.gff 2>> {log}
            barrnap --kingdom arc --threads {threads} --quiet "$f" > {output}/"$b".arc.gff 2>> {log}
            tRNAscan-SE -G -Q -q --thread {threads} -o {output}/"$b".trna.txt "$f" >> {log} 2>&1
        done
        """


checkpoint select_mags:
    input:
        bins=f"{OUT}/bins/candidates",
        checkm2=f"{OUT}/bins/checkm2/quality_report.tsv",
        gunc=f"{OUT}/bins/gunc",
        rrna_trna=f"{OUT}/bins/rrna_trna",
    output:
        mags=directory(f"{OUT}/mags/fasta"),
        summary=f"{OUT}/mags/mag_summary.tsv",
        crosswalk=f"{OUT}/mags/MAG_CROSSWALK.tsv",
    params:
        prefix=config["project_prefix"],
        min_comp=config["mag_filter"]["min_completeness"],
        max_cont=config["mag_filter"]["max_contamination"],
        gunc_required=config["mag_filter"]["require_gunc_pass"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/select_mags.py"


def mag_ids(wc=None):
    d = checkpoints.select_mags.get().output.mags
    return sorted(glob_wildcards(os.path.join(d, "{mag}.fa")).mag)
