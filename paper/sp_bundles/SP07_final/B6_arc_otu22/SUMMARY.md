# B6 — Arc_Otu22 amplicon version

**QC: PASS.** The gene aligned to itself gives 100.0%.

**Same sequence.** The server Otu22 (`Etoliko_138_2/Archaea/out/otus.fa`) and the manuscript's `Arc_Otu22_Illumina` are identical. Both are 450 bp with primers, with 0 mismatches between them. Primer-free (cutadapt linked adapter Arc519F…revcomp(Arc915R), both primers required) each is 415 bp, md5 `8bd6b191f8b4eaa35bd2d79b2126c333`.

| amplicon | md5 (primer-free) | length | identical positions | iddef0 | iddef1 | iddef2 | gene span |
|---|---|---|---|---|---|---|---|
| ms Arc_Otu22_Illumina | 8bd6b191… | 415 | 405 (405/415 = 97.6%) | 97.6% | 25.9% | 97.1% | 496–912 |
| server Otu22 | 8bd6b191… | 415 | 405 | 97.6% | 25.9% | 97.1% | 496–912 |

iddef1 is low because vsearch's edit-distance definition counts the terminal overhang of the 1,562-bp gene. It is not a measure of the amplicon region.

**Answer:** measured with vsearch, the manuscript amplicon gives **405/415 (97.6%)**, or 97.1% with the vsearch default (iddef2). It does **not** give 411/415. The 97.1% in SP03 came from this same sequence. The gap to 411/415 is aligner-dependent: in the 841–850 window (`B6_window_841_850.tsv`), vsearch aligns a stretch of the medaka-polished gene as 8 mismatches, while BLAST (SP06) resolves it as single-base indels and reaches 411/415. The 411/415 figure therefore holds only for BLAST's alignment of that window.

Files: `B6_result.tsv`, `B6_window_841_850.tsv`, `QC.txt`, `work/*.aln`, `COMMANDS.log`, `versions.tsv`, `inputs_md5.tsv`.
