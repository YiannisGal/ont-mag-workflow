# ---------------------------------------------------------------------
# Read merging, de-duplication, length filtering, statistics
# ---------------------------------------------------------------------


def unit_fastqs(sample, profile_only=False):
    u = UNITS[UNITS["sample"] == sample]
    if profile_only:
        u = u[u["profile"] == "yes"]
    return u["fastq"].tolist()


DEDUP = "seqkit rmdup -n -j {threads}" if config["reads"]["dedup_by_name"] else "seqkit seq -j {threads}"


rule merge_reads:
    """All runs of a sample, de-duplicated by read ID (used for polishing, mapping, detection)."""
    input:
        lambda wc: unit_fastqs(wc.sample),
    output:
        f"{OUT}/reads/{{sample}}.all.fastq.gz",
    log:
        f"{OUT}/logs/reads/{{sample}}.all.log",
    threads: T["small"]
    conda:
        "../envs/seqkit.yaml"
    shell:
        "(cat {input} | " + DEDUP + " -o {output}) 2> {log}"


rule merge_profile_reads:
    """Only runs flagged profile=yes (e.g. excludes amplicon runs from read-level analyses)."""
    input:
        lambda wc: unit_fastqs(wc.sample, profile_only=True),
    output:
        f"{OUT}/reads/{{sample}}.profile.fastq.gz",
    log:
        f"{OUT}/logs/reads/{{sample}}.profile.log",
    threads: T["small"]
    conda:
        "../envs/seqkit.yaml"
    shell:
        "(cat {input} | " + DEDUP + " -o {output}) 2> {log}"


rule assembly_reads:
    input:
        f"{OUT}/reads/{{sample}}.all.fastq.gz",
    output:
        f"{OUT}/reads/{{sample}}.asm.fastq.gz",
    params:
        minlen=config["reads"]["min_len_assembly"],
    threads: T["small"]
    conda:
        "../envs/seqkit.yaml"
    shell:
        "seqkit seq -j {threads} -m {params.minlen} {input} -o {output}"


rule read_stats:
    input:
        expand(f"{OUT}/reads/{{s}}.all.fastq.gz", s=ALL),
        expand(f"{OUT}/reads/{{s}}.asm.fastq.gz", s=ASM),
        expand(f"{OUT}/reads/{{s}}.profile.fastq.gz", s=PROFILE),
    output:
        f"{OUT}/qc/read_stats.tsv",
    threads: T["small"]
    conda:
        "../envs/seqkit.yaml"
    shell:
        "seqkit stats -a -T -j {threads} {input} > {output}"
