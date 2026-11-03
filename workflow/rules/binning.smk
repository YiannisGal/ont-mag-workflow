# ---------------------------------------------------------------------
# Three binners on differential coverage, reconciled with DAS_Tool
# Each binner writes <asm>/<binner>.contig2bin.tsv  (contig <tab> bin)
# ---------------------------------------------------------------------
BIN = config["binning"]


rule metabat2:
    input:
        asm=f"{OUT}/assembly/{{asm}}.contigs.fa",
        bams=binning_bams,
    output:
        depth=f"{OUT}/binning/{{asm}}/metabat2/depth.txt",
        c2b=f"{OUT}/binning/{{asm}}/metabat2.contig2bin.tsv",
    params:
        outdir=f"{OUT}/binning/{{asm}}/metabat2/bins",
        minlen=BIN["min_contig_len"],
        extra=BIN["metabat2_extra"],
    log:
        f"{OUT}/logs/binning/{{asm}}.metabat2.log",
    threads: T["medium"]
    conda:
        "../envs/metabat2.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bams} > {log} 2>&1
        metabat2 -i {input.asm} -a {output.depth} -o {params.outdir}/bin \
            -t {threads} -m {params.minlen} {params.extra} >> {log} 2>&1
        shopt -s nullglob
        for f in {params.outdir}/*.fa; do
            b=$(basename "$f" .fa | tr '.' '_')
            grep '>' "$f" | sed -e 's/^>//' -e 's/[[:space:]].*//' \
                | awk -v b="metabat2_$b" '{{print $0"\\t"b}}'
        done > {output.c2b}
        """


rule vamb:
    input:
        asm=f"{OUT}/assembly/{{asm}}.contigs.fa",
        bams=binning_bams,
    output:
        c2b=f"{OUT}/binning/{{asm}}/vamb.contig2bin.tsv",
    params:
        run=f"{OUT}/binning/{{asm}}/vamb_run",
        bamdir=f"{OUT}/mapping/{{asm}}",
        minlen=BIN["min_contig_len"],
        minfasta=BIN["min_bin_size"],
        extra=BIN["vamb_extra"],
    log:
        f"{OUT}/logs/binning/{{asm}}.vamb.log",
    threads: T["medium"]
    conda:
        "../envs/vamb.yaml"
    shell:
        """
        rm -rf {params.run}
        vamb bin default --outdir {params.run} --fasta {input.asm} --bamdir {params.bamdir} \
            -p {threads} -m {params.minlen} --minfasta {params.minfasta} {params.extra} > {log} 2>&1
        f={params.run}/vae_clusters_unsplit.tsv
        if [ -s {params.run}/vae_clusters_split.tsv ]; then f={params.run}/vae_clusters_split.tsv; fi
        awk -F'\\t' 'NR>1 {{gsub(/\\./,"_",$1); print $2"\\tvamb_"$1}}' "$f" > {output.c2b}
        """


rule semibin2:
    input:
        asm=f"{OUT}/assembly/{{asm}}.contigs.fa",
        bams=binning_bams,
    output:
        c2b=f"{OUT}/binning/{{asm}}/semibin2.contig2bin.tsv",
    params:
        outdir=f"{OUT}/binning/{{asm}}/semibin2",
        minlen=BIN["min_contig_len"],
        extra=BIN["semibin2_extra"],
    log:
        f"{OUT}/logs/binning/{{asm}}.semibin2.log",
    threads: T["medium"]
    conda:
        "../envs/semibin2.yaml"
    shell:
        """
        SemiBin2 single_easy_bin -i {input.asm} -b {input.bams} -o {params.outdir} \
            -t {threads} --min-len {params.minlen} {params.extra} > {log} 2>&1
        awk -F'\\t' 'NR>1 && $2!="-1" {{print $1"\\tsemibin2_"$2}}' {params.outdir}/contig_bins.tsv > {output.c2b}
        """


rule dastool:
    input:
        asm=f"{OUT}/assembly/{{asm}}.contigs.fa",
        metabat2=f"{OUT}/binning/{{asm}}/metabat2.contig2bin.tsv",
        vamb=f"{OUT}/binning/{{asm}}/vamb.contig2bin.tsv",
        semibin2=f"{OUT}/binning/{{asm}}/semibin2.contig2bin.tsv",
    output:
        bins=directory(f"{OUT}/binning/{{asm}}/dastool/{{asm}}_DASTool_bins"),
    params:
        prefix=f"{OUT}/binning/{{asm}}/dastool/{{asm}}",
        thr=BIN["dastool_score_threshold"],
    log:
        f"{OUT}/logs/binning/{{asm}}.dastool.log",
    threads: T["medium"]
    conda:
        "../envs/dastool.yaml"
    shell:
        """
        I=""; L=""
        for p in metabat2:{input.metabat2} vamb:{input.vamb} semibin2:{input.semibin2}; do
            l=${{p%%:*}}; f=${{p#*:}}
            if [ -s "$f" ]; then I="$I,$f"; L="$L,$l"; fi
        done
        I=${{I#,}}; L=${{L#,}}
        if [ -z "$I" ]; then
            echo "No binner produced bins for {wildcards.asm}" > {log}
            mkdir -p {output.bins}; exit 0
        fi
        DAS_Tool -i "$I" -l "$L" -c {input.asm} -o {params.prefix} --write_bins \
            -t {threads} --score_threshold {params.thr} > {log} 2>&1
        mkdir -p {output.bins}
        """
