# ---------------------------------------------------------------------
# Read-level community profiles (only runs flagged profile=yes)
# ---------------------------------------------------------------------
RP = config["read_profiling"]


rule singlem:
    input:
        f"{OUT}/reads/{{sample}}.profile.fastq.gz",
    output:
        profile=f"{OUT}/read_profiling/singlem/{{sample}}.profile.tsv",
        otus=f"{OUT}/read_profiling/singlem/{{sample}}.otu_table.tsv",
    params:
        mp=(f"--metapackage {DB['singlem_metapackage']}" if is_set(DB.get("singlem_metapackage")) else ""),
    log:
        f"{OUT}/logs/read_profiling/{{sample}}.singlem.log",
    threads: T["medium"]
    conda:
        "../envs/singlem.yaml"
    shell:
        "singlem pipe --forward {input} --taxonomic-profile {output.profile} --otu-table {output.otus} "
        "--threads {threads} {params.mp} > {log} 2>&1"


rule singlem_phylum:
    """Phylum relative abundance from filled coverage (taxon + descendants)."""
    input:
        f"{OUT}/read_profiling/singlem/{{sample}}.profile.tsv",
    output:
        f"{OUT}/read_profiling/singlem/{{sample}}.phylum.tsv",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/singlem_phylum.py"


rule kraken2:
    input:
        f"{OUT}/reads/{{sample}}.profile.fastq.gz",
    output:
        report=f"{OUT}/read_profiling/kraken2/{{sample}}.report.txt",
        bracken=f"{OUT}/read_profiling/kraken2/{{sample}}.bracken.tsv",
    params:
        db=DB.get("kraken2", ""),
        rl=RP["bracken_read_len"],
        lvl=RP["bracken_level"],
    log:
        f"{OUT}/logs/read_profiling/{{sample}}.kraken2.log",
    threads: T["medium"]
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        kraken2 --db {params.db} --threads {threads} --report {output.report} \
            --output /dev/null --gzip-compressed {input} > {log} 2>&1
        bracken -d {params.db} -i {output.report} -o {output.bracken} \
            -r {params.rl} -l {params.lvl} >> {log} 2>&1
        """


rule kaiju:
    input:
        f"{OUT}/reads/{{sample}}.profile.fastq.gz",
    output:
        raw=f"{OUT}/read_profiling/kaiju/{{sample}}.kaiju.out",
        table=f"{OUT}/read_profiling/kaiju/{{sample}}.table.tsv",
    params:
        nodes=DB.get("kaiju_nodes", ""),
        names=DB.get("kaiju_names", ""),
        fmi=DB.get("kaiju_fmi", ""),
        rank=RP["kaiju_rank"],
    log:
        f"{OUT}/logs/read_profiling/{{sample}}.kaiju.log",
    threads: T["medium"]
    conda:
        "../envs/kaiju.yaml"
    shell:
        """
        kaiju -t {params.nodes} -f {params.fmi} -i {input} -o {output.raw} -z {threads} > {log} 2>&1
        kaiju2table -t {params.nodes} -n {params.names} -r {params.rank} \
            -o {output.table} {output.raw} >> {log} 2>&1
        """


# ---------------------------------------------------------------------
# Competitive gene-panel read search (e.g. mcrA vs mtmB/mtbB/mttB/mtaB)
# ---------------------------------------------------------------------
GP = config["gene_panel"]


rule gene_panel_db:
    input:
        GP["proteins"] if is_set(GP.get("proteins")) else [],
    output:
        f"{OUT}/gene_panel/panel.dmnd",
    conda:
        "../envs/diamond.yaml"
    shell:
        "diamond makedb --in {input} -d {output}"


rule gene_panel_search:
    input:
        db=f"{OUT}/gene_panel/panel.dmnd",
        reads=f"{OUT}/reads/{{sample}}.profile.fastq.gz",
    output:
        f"{OUT}/gene_panel/{{sample}}.diamond.tsv",
    log:
        f"{OUT}/logs/gene_panel/{{sample}}.log",
    threads: T["medium"]
    conda:
        "../envs/diamond.yaml"
    shell:
        "diamond blastx -q {input.reads} -d {input.db} -o {output} -p {threads} --long-reads "
        "--outfmt 6 qseqid sseqid pident length evalue bitscore > {log} 2>&1"


rule gene_panel_counts:
    input:
        hits=f"{OUT}/gene_panel/{{sample}}.diamond.tsv",
        labels=GP["labels"] if is_set(GP.get("labels")) else [],
    output:
        f"{OUT}/gene_panel/{{sample}}.counts.tsv",
    params:
        min_bits=GP["min_bitscore"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/gene_panel_counts.py"
