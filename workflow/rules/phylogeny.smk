# ---------------------------------------------------------------------
# Phylogenomics: GTDB-Tk bac120 / ar53 marker alignment of MAGs (+ your
# reference genomes, e.g. isolates and type strains) -> IQ-TREE
# (ModelFinder, 1000 SH-aLRT + 1000 ultrafast bootstrap)
# ---------------------------------------------------------------------


rule tree_genomes:
    input:
        f"{OUT}/mags/fasta",
    output:
        directory(f"{OUT}/phylogeny/genomes"),
    params:
        refs=REFS,
    shell:
        """
        mkdir -p {output}
        shopt -s nullglob
        for f in {input}/*.fa; do ln -sf "$(realpath "$f")" {output}/; done
        if [ -n "{params.refs}" ]; then
            for f in {params.refs}/*.fa {params.refs}/*.fna {params.refs}/*.fasta; do
                n=$(basename "$f"); n=${{n%.*}}
                ln -sf "$(realpath "$f")" {output}/"$n".fa
            done
        fi
        """


rule gtdbtk_align:
    input:
        f"{OUT}/phylogeny/genomes",
    output:
        f"{OUT}/phylogeny/gtdbtk_align.done",
    params:
        identify=f"{OUT}/phylogeny/identify",
        align=f"{OUT}/phylogeny/align",
        db=DB["gtdbtk"],
    log:
        f"{OUT}/logs/phylogeny/gtdbtk_align.log",
    threads: T["big"]
    resources:
        mem_mb=MEM["gtdbtk"],
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        """
        export GTDBTK_DATA_PATH={params.db}
        gtdbtk identify --genome_dir {input} -x fa --out_dir {params.identify} --cpus {threads} > {log} 2>&1
        gtdbtk align --identify_dir {params.identify} --out_dir {params.align} \
            --skip_gtdb_refs --cpus {threads} >> {log} 2>&1
        touch {output}
        """


rule iqtree:
    input:
        f"{OUT}/phylogeny/gtdbtk_align.done",
    output:
        f"{OUT}/phylogeny/iqtree.done",
    params:
        align=f"{OUT}/phylogeny/align/align",
        outdir=f"{OUT}/phylogeny/trees",
    log:
        f"{OUT}/logs/phylogeny/iqtree.log",
    threads: T["medium"]
    conda:
        "../envs/iqtree.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        for dom in bac120 ar53; do
            msa={params.align}/gtdbtk.$dom.user_msa.fasta.gz
            [ -s "$msa" ] || continue
            n=$(zcat "$msa" | grep -c '>' || true)
            if [ "$n" -lt 4 ]; then echo "$dom: only $n genomes, no tree" >> {log}; continue; fi
            zcat "$msa" > {params.outdir}/$dom.msa.fasta
            iqtree -s {params.outdir}/$dom.msa.fasta -m MFP -B 1000 -alrt 1000 \
                -T {threads} --prefix {params.outdir}/$dom -redo >> {log} 2>&1
        done
        touch {output}
        """
