# SP06: which Arc_Otu22 amplicon gives 97.1%?

1. **Is it the same sequence?** Yes. The server's Otu22 (`Etoliko_138_2/Archaea/out/otus.fa`) and the manuscript's `Arc_Otu22_Illumina` are identical: 450 bp each, 100% identity, 0 mismatches, 0 gaps. The manuscript file was taken from the server copy in `98_arc22_longonly/` (md5 a6235994…). The gene file is md5 2d17845f…, the same as the laptop copy.

2. **What does each amplicon give against the gene?** Both give identical numbers, because they are the same sequence.

| method | identity | detail |
|---|---|---|
| vsearch iddef 0 | 97.6% | 452 aligned columns, 11 mismatches, 2 gaps |
| vsearch iddef 1 | no hit at `--id 0.5` | edit-distance identity is computed over the whole 1,562-bp gene |
| vsearch iddef 2 (default) | 97.1% | same alignment as iddef 0 |
| **BLASTN, primers removed (415 bp)** | **411/415** | 418 aligned columns, 3 mismatches, 4 gaps (gene positions 496–912) |

3. **Verdict: same sequence.** The difference is not the OTU table and not mainly the identity definition. It comes from two things:
   - **Alignment method.** vsearch's global aligner scores the divergent stretch at gene positions ~841–850 as a run of mismatches. BLAST places single-base gaps there instead, which matches the three extra insertions in the gene.
   - **Primer handling.** vsearch was run on the 450-bp amplicon with primers still attached; BLAST on the 415-bp primer-free template.

   **The thesis keeps 411/415 (99.0%), with the note that the full-length gene carries three single-base insertions.**

Commands and versions: `COMMANDS.log` (vsearch 2.32.0, BLAST+ 2.12.0). Results: `SP06_RESULT.tsv`.
