# ---------------------------------------------------------------------
# Per-MAG annotation (runs after the select_mags checkpoint)
# ---------------------------------------------------------------------


rule prodigal:
    input:
        f"{OUT}/mags/fasta/{{mag}}.fa",
    output:
        faa=f"{OUT}/annotation/prodigal/{{mag}}.faa",
        fna=f"{OUT}/annotation/prodigal/{{mag}}.fna",
        gff=f"{OUT}/annotation/prodigal/{{mag}}.gff",
    conda:
        "../envs/prodigal.yaml"
    shell:
        "prodigal -i {input} -a {output.faa} -d {output.fna} -f gff -o {output.gff} -p single -q"


rule kofamscan:
    """KEGG orthologs with the curated, profile-specific thresholds (only '*' rows kept).
    A flat E-value cut-off grossly inflates hits, especially for transporters."""
    input:
        f"{OUT}/annotation/prodigal/{{mag}}.faa",
    output:
        full=f"{OUT}/annotation/kofam/{{mag}}.detail.tsv",
        sig=f"{OUT}/annotation/kofam/{{mag}}.significant.tsv",
    params:
        profiles=DB["kofam_profiles"],
        ko_list=DB["kofam_ko_list"],
        tmp=f"{OUT}/annotation/kofam/tmp_{{mag}}",
    log:
        f"{OUT}/logs/annotation/{{mag}}.kofam.log",
    threads: T["medium"]
    conda:
        "../envs/kofamscan.yaml"
    shell:
        """
        exec_annotation -o {output.full} --profile {params.profiles} --ko-list {params.ko_list} \
            --cpu {threads} --format detail-tsv --tmp-dir {params.tmp} {input} > {log} 2>&1
        rm -rf {params.tmp}
        (grep '^\\*' {output.full} || true) > {output.sig}
        """


rule dbcan:
    input:
        f"{OUT}/annotation/prodigal/{{mag}}.faa",
    output:
        f"{OUT}/annotation/dbcan/{{mag}}/overview.txt",
    params:
        outdir=f"{OUT}/annotation/dbcan/{{mag}}",
        db=DB["dbcan"],
    log:
        f"{OUT}/logs/annotation/{{mag}}.dbcan.log",
    threads: T["small"]
    conda:
        "../envs/dbcan.yaml"
    shell:
        "run_dbcan {input} protein --out_dir {params.outdir} --db_dir {params.db} "
        "--dia_cpu {threads} --hmm_cpu {threads} --tf_cpu {threads} > {log} 2>&1"


rule genomad:
    input:
        f"{OUT}/mags/fasta/{{mag}}.fa",
    output:
        directory(f"{OUT}/annotation/genomad/{{mag}}"),
    params:
        db=DB["genomad"],
    log:
        f"{OUT}/logs/annotation/{{mag}}.genomad.log",
    threads: T["small"]
    conda:
        "../envs/genomad.yaml"
    shell:
        "genomad end-to-end --cleanup --threads {threads} {input} {output} {params.db} > {log} 2>&1"


rule antismash:
    input:
        f"{OUT}/mags/fasta/{{mag}}.fa",
    output:
        directory(f"{OUT}/annotation/antismash/{{mag}}"),
    params:
        db=(f"--databases {DB['antismash']}" if is_set(DB.get("antismash")) else ""),
    log:
        f"{OUT}/logs/annotation/{{mag}}.antismash.log",
    threads: T["small"]
    conda:
        "../envs/antismash.yaml"
    shell:
        "antismash --genefinding-tool prodigal-m --cb-knownclusters --cpus {threads} {params.db} "
        "--output-dir {output} {input} > {log} 2>&1"


rule bakta:
    input:
        f"{OUT}/mags/fasta/{{mag}}.fa",
    output:
        directory(f"{OUT}/annotation/bakta/{{mag}}"),
    params:
        db=DB.get("bakta", ""),
    log:
        f"{OUT}/logs/annotation/{{mag}}.bakta.log",
    threads: T["small"]
    conda:
        "../envs/bakta.yaml"
    shell:
        "bakta --db {params.db} --output {output} --prefix {wildcards.mag} "
        "--locus-tag {wildcards.mag} --threads {threads} --meta {input} > {log} 2>&1"


# ---------------- aggregation ----------------------------------------


rule ko_matrix:
    input:
        lambda wc: expand(f"{OUT}/annotation/kofam/{{mag}}.significant.tsv", mag=mag_ids()),
    output:
        f"{OUT}/annotation/summary/ko_matrix.tsv",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/ko_matrix.py"


rule cazy_matrix:
    input:
        lambda wc: expand(f"{OUT}/annotation/dbcan/{{mag}}/overview.txt", mag=mag_ids()),
    output:
        f"{OUT}/annotation/summary/cazy_matrix.tsv",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/cazy_matrix.py"


rule genomad_done:
    input:
        lambda wc: expand(f"{OUT}/annotation/genomad/{{mag}}", mag=mag_ids()),
    output:
        f"{OUT}/annotation/summary/genomad_done.txt",
    shell:
        "ls -d {input} > {output}"


rule antismash_done:
    input:
        lambda wc: expand(f"{OUT}/annotation/antismash/{{mag}}", mag=mag_ids()),
    output:
        f"{OUT}/annotation/summary/antismash_done.txt",
    shell:
        "ls -d {input} > {output}"


rule bakta_done:
    input:
        lambda wc: expand(f"{OUT}/annotation/bakta/{{mag}}", mag=mag_ids()),
    output:
        f"{OUT}/annotation/summary/bakta_done.txt",
    shell:
        "ls -d {input} > {output}"


# ---------------- metabolic models (optional) ------------------------


rule gapseq:
    input:
        f"{OUT}/mags/fasta/{{mag}}.fa",
    output:
        f"{OUT}/metabolic_models/gapseq/{{mag}}/{{mag}}.xml",
    params:
        outdir=f"{OUT}/metabolic_models/gapseq/{{mag}}",
    log:
        f"{OUT}/logs/metabolic_models/{{mag}}.gapseq.log",
    threads: T["small"]
    conda:
        "../envs/gapseq.yaml"
    shell:
        """
        fa=$(realpath {input}); lg=$(realpath -m {log})
        cd {params.outdir}
        gapseq doall "$fa" -K {threads} > "$lg" 2>&1
        """


rule smetana:
    input:
        lambda wc: expand(f"{OUT}/metabolic_models/gapseq/{{mag}}/{{mag}}.xml", mag=mag_ids()),
    output:
        f"{OUT}/metabolic_models/smetana/smetana_detailed.tsv",
    params:
        prefix=f"{OUT}/metabolic_models/smetana/smetana",
    log:
        f"{OUT}/logs/metabolic_models/smetana.log",
    conda:
        "../envs/smetana.yaml"
    shell:
        "smetana {input} -d --flavor fbc2 -o {params.prefix} > {log} 2>&1"
