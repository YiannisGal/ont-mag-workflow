# ---------------------------------------------------------------------
# Every 'assemble' read set is mapped to every assembly
# (self + cross mapping -> differential coverage for binning)
# ---------------------------------------------------------------------


rule map_to_assembly:
    input:
        asm=f"{OUT}/assembly/{{asm}}.contigs.fa",
        reads=f"{OUT}/reads/{{reads}}.all.fastq.gz",
    output:
        bam=f"{OUT}/mapping/{{asm}}/{{reads}}.bam",
        bai=f"{OUT}/mapping/{{asm}}/{{reads}}.bam.bai",
        flagstat=f"{OUT}/mapping/{{asm}}/flagstat/{{reads}}.flagstat.txt",
    log:
        f"{OUT}/logs/mapping/{{asm}}.{{reads}}.log",
    threads: T["medium"]
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        (minimap2 -ax map-ont -t {threads} {input.asm} {input.reads} \
            | samtools sort -@ {threads} -o {output.bam} -) 2> {log}
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        """


def binning_bams(wc):
    return expand(f"{OUT}/mapping/{wc.asm}/{{r}}.bam", r=ASM)
