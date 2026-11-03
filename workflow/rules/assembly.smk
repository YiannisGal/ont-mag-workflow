# ---------------------------------------------------------------------
# Flye metagenome assembly (reads >= min length) + Medaka polishing (all reads)
# ---------------------------------------------------------------------


rule flye:
    input:
        f"{OUT}/reads/{{sample}}.asm.fastq.gz",
    output:
        f"{OUT}/assembly/{{sample}}/flye/assembly.fasta",
    params:
        outdir=f"{OUT}/assembly/{{sample}}/flye",
        rt=config["assembly"]["flye_read_type"],
        extra=config["assembly"]["flye_extra"],
    log:
        f"{OUT}/logs/assembly/{{sample}}.flye.log",
    threads: T["big"]
    resources:
        mem_mb=MEM["assembly"],
    conda:
        "../envs/flye.yaml"
    shell:
        "flye {params.rt} {input} --out-dir {params.outdir} --threads {threads} {params.extra} > {log} 2>&1"


rule medaka:
    input:
        reads=f"{OUT}/reads/{{sample}}.all.fastq.gz",
        draft=f"{OUT}/assembly/{{sample}}/flye/assembly.fasta",
    output:
        f"{OUT}/assembly/{{sample}}/medaka/consensus.fasta",
    params:
        outdir=f"{OUT}/assembly/{{sample}}/medaka",
        model=config["assembly"]["medaka_model"],
    log:
        f"{OUT}/logs/assembly/{{sample}}.medaka.log",
    threads: T["medium"]
    conda:
        "../envs/medaka.yaml"
    shell:
        "medaka_consensus -i {input.reads} -d {input.draft} -o {params.outdir} "
        "-t {threads} -m {params.model} -f > {log} 2>&1"


def final_draft(wc):
    if config["assembly"].get("polish", True):
        return f"{OUT}/assembly/{wc.sample}/medaka/consensus.fasta"
    return f"{OUT}/assembly/{wc.sample}/flye/assembly.fasta"


rule rename_contigs:
    """Prefix contigs with the sample name (SED_contig_1 ...) so assemblies never collide."""
    input:
        final_draft,
    output:
        f"{OUT}/assembly/{{sample}}.contigs.fa",
    conda:
        "../envs/seqkit.yaml"
    shell:
        "seqkit replace -p '^(\\S+).*$' -r '{wildcards.sample}_${{1}}' {input} > {output}"


rule assembly_stats:
    input:
        expand(f"{OUT}/assembly/{{s}}.contigs.fa", s=ASM),
    output:
        f"{OUT}/qc/assembly_stats.tsv",
    conda:
        "../envs/seqkit.yaml"
    shell:
        "seqkit stats -a -T {input} > {output}"
