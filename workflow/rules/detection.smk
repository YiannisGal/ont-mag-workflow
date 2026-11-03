# ---------------------------------------------------------------------
# Competitive recruitment of every read set (incl. detect_only samples,
# e.g. public metagenomes) to all MAGs (+ optional background genomes).
# Outputs: relative abundance / mean depth / breadth, a breadth sweep over
# read-identity thresholds, and a same-species / related / not-detected call.
# ---------------------------------------------------------------------
DET = config["detection"]
BG = DET.get("background_genomes", "")


rule detection_reference:
    input:
        f"{OUT}/mags/fasta",
    output:
        f"{OUT}/detection/reference.fa",
    params:
        bg=BG,
    shell:
        """
        shopt -s nullglob
        files=({input}/*.fa)
        if [ -n "{params.bg}" ]; then files+=({params.bg}/*.fa {params.bg}/*.fna {params.bg}/*.fasta); fi
        for f in "${{files[@]}}"; do
            n=$(basename "$f"); n=${{n%.*}}
            awk -v g="$n" '/^>/ {{sub(/^>/, ">" g "~"); sub(/[ \\t].*/, "")}} 1' "$f"
        done > {output}
        """


rule detection_map:
    input:
        ref=f"{OUT}/detection/reference.fa",
        reads=f"{OUT}/reads/{{reads}}.all.fastq.gz",
    output:
        bam=f"{OUT}/detection/bam/{{reads}}.bam",
        bai=f"{OUT}/detection/bam/{{reads}}.bam.bai",
    log:
        f"{OUT}/logs/detection/{{reads}}.map.log",
    threads: T["medium"]
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        (minimap2 -ax map-ont --secondary=no -t {threads} {input.ref} {input.reads} \
            | samtools sort -@ {threads} -o {output.bam} -) 2> {log}
        samtools index {output.bam}
        """


rule coverm_base:
    input:
        f"{OUT}/detection/bam/{{reads}}.bam",
    output:
        f"{OUT}/detection/coverm/{{reads}}.base.tsv",
    threads: T["small"]
    conda:
        "../envs/coverm.yaml"
    shell:
        "coverm genome -b {input} -s '~' -t {threads} --min-covered-fraction 0 "
        "-m relative_abundance mean covered_fraction > {output}"


rule coverm_identity:
    input:
        f"{OUT}/detection/bam/{{reads}}.bam",
    output:
        f"{OUT}/detection/coverm/{{reads}}.id{{ident}}.tsv",
    params:
        aln=DET["min_read_aligned_percent"],
    threads: T["small"]
    conda:
        "../envs/coverm.yaml"
    shell:
        "coverm genome -b {input} -s '~' -t {threads} --min-covered-fraction 0 "
        "--min-read-percent-identity {wildcards.ident} --min-read-aligned-percent {params.aln} "
        "-m covered_fraction mean > {output}"


rule detection_calls:
    input:
        bam=f"{OUT}/detection/bam/{{reads}}.bam",
        base=f"{OUT}/detection/coverm/{{reads}}.base.tsv",
        ids=expand(f"{OUT}/detection/coverm/{{{{reads}}}}.id{{i}}.tsv", i=DET["identities"]),
    output:
        f"{OUT}/detection/{{reads}}.calls.tsv",
    params:
        identities=DET["identities"],
        min_aln=DET["min_read_aligned_percent"],
        mode_thr=DET["same_species_identity_mode"],
        same_breadth=DET["same_species_breadth"],
        related_breadth=DET["related_breadth"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/detection_calls.py"


rule detection_merge:
    input:
        expand(f"{OUT}/detection/{{reads}}.calls.tsv", reads=ALL),
    output:
        f"{OUT}/detection/coverm_all.tsv",
    run:
        dfs = []
        for f in input:
            d = pd.read_csv(f, sep="\t")
            d.insert(0, "read_set", os.path.basename(f).replace(".calls.tsv", ""))
            dfs.append(d)
        pd.concat(dfs).to_csv(output[0], sep="\t", index=False)
