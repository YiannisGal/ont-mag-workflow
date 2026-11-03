#!/usr/bin/env bash
# 56_dsr_placement.sh -- identify the two unbinned dsrAB organisms, properly.
#
#   bash scripts/run.sh scripts/56_dsr_placement_XXXXXXX.sh
#
# WHY THIS EXISTS. Step 55 returned a strikingly clean read-level answer --
# LIQ_contig_12880 reads called Desulfosporosinus 28/33, LIQ_contig_19299 reads
# called Peptococcaceae across several genera -- AND A VOID CONTROL. The control
# contig came back with zero reads because I built the mapping reference from
# the two dsr contigs ALONE, so no read could ever map to a control contig. The
# control did not fail; it was never run. 55 section E is therefore unquotable
# as it stands.
#
# WHAT THIS DOES.
#   A  Rebuilds the mapping reference with the control contigs IN IT, and runs
#      three controls chosen to be able to fail:
#        - a Bacillota MAG with a known GTDB genus
#        - an archaeal MAG (must NOT come back as a sulfate reducer)
#        - LIQ_das_metabat_022, the Desulfohalotomaculaceae MAG. This is the
#          sharp one: if Kaiju calls ITS reads Peptococcaceae too, then
#          "Peptococcaceae" is what Kaiju says about everything in this family
#          and the dsr result means less than it looks.
#   B  Fetches DsrA and DsrB reference proteins from NCBI, the same way step 15
#      fetched its MCR references, spanning the reductive lineages AND the
#      reverse-type (sulfur-oxidiser) dsr that must be excluded.
#   C  blastp: nearest named relatives of each query protein, with identity.
#   D  mafft + FastTree, and the sister group of each query read off the tree.
#
# Read-only except 56_dsr_placement/.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/56_dsr_placement"; mkdir -p "$OUT"
export TMPDIR="${TMPDIR:-$WS2/tmp}"; mkdir -p "$TMPDIR"
exec > >(tee "$OUT/56.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }; no(){ printf '  [--]   %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
have(){ command -v "$1" >/dev/null 2>&1; }

echo "56 -- dsrAB placement"
echo "date : $(date -Is)"
echo "out  : $OUT"

ASM="$WS2/04_polish/LIQ_asm_final.fasta"
DSRFAA="$WS2/49_dsr/dsr_contig_proteins.faa"
KJ="$WS2/35_kaiju/LIQ.kaiju"
NAMES="$HOME/databases/kraken2_std16/names.dmp"
BINS="$WS2/06_binqc/union_bins"
GTDB_B="$WS2/08_taxonomy/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
GTDB_A="$WS2/08_taxonomy/gtdbtk/classify/gtdbtk.ar53.summary.tsv"
DSR_CONTIGS="LIQ_contig_12880 LIQ_contig_19299"

for f in "$ASM" "$DSRFAA" "$KJ"; do
  [ -s "$f" ] || { bad "missing: $f"; echo "STEP_STATUS: FAILED"; exit 1; }
done
[ -s "$NAMES" ] || no "names.dmp missing -- taxids will print as numbers"

# ================================================================ A
say "A. THE CONTROL, THIS TIME ACTUALLY IN THE REFERENCE"

pick_contigs(){   # $1 = bin name, $2 = how many
  local B="$1" N="$2" FA
  FA="$(find "$BINS" -maxdepth 1 -name "${B}.fa*" | head -1)"
  [ -z "$FA" ] && return
  seqkit fx2tab -nl "$FA" 2>/dev/null | sort -k2,2nr | head -"$N" \
    | awk '{gsub(/.*~/,"",$1); print $1}'
}

CTRL_BACT="LIQ_das_metabat_001"
CTRL_SRB="LIQ_das_metabat_022"
CTRL_ARC="LIQ_das_vamb_017"

sub "what GTDB-Tk already says about each control"
for B in "$CTRL_BACT" "$CTRL_SRB" "$CTRL_ARC"; do
  T="$(awk -F'\t' -v b="$B" 'NR>1 && $1==b{print $2}' "$GTDB_B" "$GTDB_A" 2>/dev/null | head -1)"
  printf '  %-24s %s\n' "$B" "${T:-<not in either summary>}"
done

CTRL_CONTIGS=""
for B in "$CTRL_BACT" "$CTRL_SRB" "$CTRL_ARC"; do
  C="$(pick_contigs "$B" 2 | tr '\n' ' ')"
  printf '  %-24s contigs: %s\n' "$B" "${C:-NONE FOUND}"
  CTRL_CONTIGS="$CTRL_CONTIGS $C"
done

ALL_CONTIGS="$DSR_CONTIGS $CTRL_CONTIGS"
samtools faidx "$ASM" $ALL_CONTIGS > "$OUT/ref_contigs.fasta" 2>"$OUT/faidx.err"
NC=$(grep -c '^>' "$OUT/ref_contigs.fasta" 2>/dev/null || echo 0)
if [ "$NC" -lt 3 ]; then
  bad "only $NC contigs extracted -- see faidx.err"; echo "STEP_STATUS: FAILED"; exit 1
fi
ok "$NC contigs in the mapping reference (2 dsr + $((NC-2)) control)"

READS="$WS2/00_input/liquid_reads.fastq.gz"
if [ ! -s "$OUT/ctrl_reads.bam" ]; then
  echo "  mapping $(basename "$READS") -> $NC contigs ..."
  minimap2 -x map-ont -a --secondary=no -t 8 "$OUT/ref_contigs.fasta" "$READS" 2>"$OUT/minimap2.err" \
    | samtools sort -@ 4 -o "$OUT/ctrl_reads.bam" - 2>/dev/null
  samtools index "$OUT/ctrl_reads.bam" 2>/dev/null
fi
[ -s "$OUT/ctrl_reads.bam" ] || { bad "mapping produced nothing"; echo "STEP_STATUS: FAILED"; exit 1; }
ok "${OUT#$WS2/}/ctrl_reads.bam"

lookup(){   # $1 contig  $2 label
  local C="$1" L="$2" IDS="$TMPDIR/.ids56" N M
  samtools view "$OUT/ctrl_reads.bam" "$C" 2>/dev/null | awk '{print $1}' | sort -u > "$IDS"
  N=$(wc -l < "$IDS")
  local LEN; LEN=$(awk -v c="$C" '$1==c{print $2}' "$ASM.fai")
  printf '\n  %-12s %-20s %5d reads on %6s bp\n' "$L" "$C" "$N" "${LEN:-?}"
  [ "$N" -eq 0 ] && { echo "      (none)"; return; }
  awk -F'\t' -v ids="$IDS" 'BEGIN{while((getline l<ids)>0) k[l]=1} ($2 in k){print}' "$KJ" > "$TMPDIR/.tax56"
  M=$(awk -F'\t' '$1=="C"' "$TMPDIR/.tax56" | wc -l)
  printf '      classified %d / %d (%.0f%%)\n' "$M" "$N" "$(awk -v a=$M -v b=$N 'BEGIN{print (b?100*a/b:0)}')"
  awk -F'\t' '$1=="C"{print $3}' "$TMPDIR/.tax56" | sort | uniq -c | sort -rn | head -8 \
  | while read -r cnt tid; do
      nm=""
      [ -s "$NAMES" ] && nm=$(awk -F'\t\\|\t' -v t="$tid" '$1==t && $4 ~ /scientific name/{print $2; exit}' "$NAMES")
      printf '      %5d  %-9s %s\n' "$cnt" "$tid" "${nm:-<taxid>}"
    done
}

sub "CONTROLS -- each must return its own known lineage"
i=0
for B in "$CTRL_BACT" "$CTRL_SRB" "$CTRL_ARC"; do
  i=$((i+1))
  for C in $(pick_contigs "$B" 2); do lookup "$C" "CTRL$i"; done
done

sub "THE QUESTION"
for C in $DSR_CONTIGS; do lookup "$C" "dsr"; done

# ================================================================ B
say "B. DsrA / DsrB REFERENCE PROTEINS FROM NCBI"

EUT="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
fetch(){  # $1 outfile  $2 retmax  $3 term
  local out="$1" n="$2" term="$3" enc ids
  enc=$(printf '%s' "$term" | sed 's/ /+/g; s/"/%22/g; s/\[/%5B/g; s/\]/%5D/g')
  ids=$(wget -qO- "$EUT/esearch.fcgi?db=protein&term=${enc}&retmax=${n}" 2>/dev/null \
        | grep -oP '(?<=<Id>)[0-9]+(?=</Id>)' | tr '\n' ',' | sed 's/,$//')
  if [ -z "$ids" ]; then printf '      %-58s 0\n' "$term"; return 1; fi
  wget -qO- "$EUT/efetch.fcgi?db=protein&id=${ids}&rettype=fasta&retmode=text" 2>/dev/null >> "$out"
  printf '      %-58s %s\n' "$term" "$(printf '%s' "$ids" | tr ',' '\n' | wc -l)"
  sleep 0.4
}

if ! wget -q --spider --timeout=15 "$EUT/einfo.fcgi" 2>/dev/null; then
  bad "NCBI is not reachable from this machine. Sections B-D cannot run."
  echo "         Section A stands on its own; see the verdict at the end."
  NET=0
else
  ok "NCBI reachable"
  NET=1
fi

if [ "$NET" = 1 ] && [ ! -s "$OUT/refs_raw.faa" ]; then
  : > "$OUT/refs_raw.faa"
  for SUB in "alpha" "beta"; do
    sub "DsrA/B $SUB subunit references"
    fetch "$OUT/refs_raw.faa" 80 "\"dissimilatory sulfite reductase $SUB subunit\"[Protein Name] AND refseq[filter]"
    for ORG in Desulfosporosinus Desulfotomaculum Pelotomaculum Desulfovibrio \
               Desulfobacterium Archaeoglobus Thermodesulfovibrio Desulfurivibrio \
               Allochromatium Thiobacillus; do
      fetch "$OUT/refs_raw.faa" 6 "\"dissimilatory sulfite reductase $SUB subunit\"[Protein Name] AND \"$ORG\"[Organism]"
    done
  done
fi

if [ -s "$OUT/refs_raw.faa" ]; then
  # rename to REF|Genus_species|accession, drop duplicates by accession
  awk '
    /^>/ { acc=substr($1,2)
           org=""
           if (match($0,/\[[^]]*\]$/)) org=substr($0,RSTART+1,RLENGTH-2)
           gsub(/[^A-Za-z0-9]+/,"_",org)
           if (org=="") org="unknown"
           if (acc in seen) { skip=1; next }
           seen[acc]=1; skip=0
           print ">REF|" org "|" acc
           next }
    { if(!skip) print }' "$OUT/refs_raw.faa" > "$OUT/refs.faa"
  NR_=$(grep -c '^>' "$OUT/refs.faa")
  ok "$NR_ unique reference proteins"
  echo "  genera represented:"
  grep '^>' "$OUT/refs.faa" | awk -F'|' '{split($2,a,"_"); print a[1]}' | sort | uniq -c | sort -rn | head -20 | sed 's/^/    /'
else
  no "no references retrieved"
fi

# ================================================================ C
say "C. NEAREST NAMED RELATIVES OF EACH QUERY PROTEIN"

QIDS="LIQ_contig_12880_3 LIQ_contig_12880_2 LIQ_contig_19299_14 LIQ_contig_19299_15 LIQ_contig_19299_13"
seqkit grep -n -f <(printf '%s\n' $QIDS) "$DSRFAA" 2>/dev/null \
  | seqkit replace -p '\s.*$' -r '' 2>/dev/null \
  | sed 's/^>/>QUERY|/' > "$OUT/queries.faa"
NQ=$(grep -c '^>' "$OUT/queries.faa" 2>/dev/null || echo 0)
if [ "$NQ" -eq 0 ]; then
  # fallback: plain awk extraction
  awk -v ids="$QIDS" 'BEGIN{n=split(ids,a," "); for(i=1;i<=n;i++) want[a[i]]=1}
    /^>/{ id=substr($1,2); p=(id in want); if(p) print ">QUERY|" id; next }
    { if(p) print }' "$DSRFAA" > "$OUT/queries.faa"
  NQ=$(grep -c '^>' "$OUT/queries.faa")
fi
ok "$NQ query proteins"
grep '^>' "$OUT/queries.faa" | sed 's/^/    /'

echo
echo "  expected from step 49: _12880_3 dsrA, _12880_2 dsrB, _19299_14 dsrA,"
echo "  _19299_13 dsrB. _19299_15 is included because it sits immediately"
echo "  upstream of _19299_14 and overlaps it by 4 bp -- if dsrA here is split"
echo "  the way MAG_017's mcrA is, this is the other half."

if [ -s "$OUT/refs.faa" ] && [ "$NQ" -gt 0 ]; then
  makeblastdb -in "$OUT/refs.faa" -dbtype prot -out "$OUT/refdb" >/dev/null 2>&1
  blastp -query "$OUT/queries.faa" -db "$OUT/refdb" -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' \
         -max_target_seqs 6 -num_threads 8 > "$OUT/blastp.tsv" 2>/dev/null
  if [ -s "$OUT/blastp.tsv" ]; then
    awk -F'\t' '
      { if($1!=q){ q=$1; printf "\n  %s\n", q
                   printf "      %-42s %6s %7s %6s %6s %10s\n","nearest named relative","pid%","alnlen","qlen","slen","evalue" }
        split($2,s,"|"); printf "      %-42s %6.1f %7d %6d %6d %10s\n", s[2]" ("s[3]")", $3, $4, $5, $6, $7 }' "$OUT/blastp.tsv"
  else
    no "blastp returned nothing"
  fi
else
  no "skipped -- need both references and queries"
fi

# ================================================================ D
say "D. PHYLOGENETIC PLACEMENT"

if [ -s "$OUT/refs.faa" ] && [ "$NQ" -gt 0 ]; then
  cat "$OUT/queries.faa" "$OUT/refs.faa" > "$OUT/all.faa"
  mafft --auto --quiet --thread 8 "$OUT/all.faa" > "$OUT/all.aln" 2>"$OUT/mafft.err"
  if [ -s "$OUT/all.aln" ]; then
    ok "aligned $(grep -c '^>' "$OUT/all.aln") sequences"
    FastTree -lg -quiet "$OUT/all.aln" > "$OUT/dsr.tree" 2>"$OUT/fasttree.err"
    if [ -s "$OUT/dsr.tree" ]; then
      ok "tree: ${OUT#$WS2/}/dsr.tree"
      python3 - "$OUT/dsr.tree" <<'PY'
import re, sys
s = open(sys.argv[1]).read().strip().rstrip(';')
pos = [0]
def node():
    ch = []
    if s[pos[0]] == '(':
        pos[0] += 1
        while True:
            ch.append(node())
            if s[pos[0]] == ',': pos[0] += 1
            else: break
        pos[0] += 1
    m = re.compile(r"[^,():;]*").match(s, pos[0]); lab = m.group(); pos[0] = m.end()
    if pos[0] < len(s) and s[pos[0]] == ':':
        pos[0] += 1
        m = re.compile(r"[-0-9.eE+]*").match(s, pos[0]); pos[0] = m.end()
    return {'ch': ch, 'lab': lab}
root = node()
def leaves(n): return [n['lab']] if not n['ch'] else [l for c in n['ch'] for l in leaves(c)]
allq = [l for l in leaves(root) if l.startswith('QUERY')]
def clades(target):
    path = []
    def go(n):
        if not n['ch']: return n['lab'] == target
        if any(go(c) for c in n['ch']): path.append(n); return True
        return False
    go(root)
    return [[x for x in leaves(n) if x != target] for n in path]
for q in allq:
    print("\n  %s" % q)
    shown = 0
    for cl in clades(q):
        refs = [c for c in cl if c.startswith('REF')]
        if not refs: continue
        gen = {}
        for r in refs:
            g = r.split('|')[1].split('_')[0]
            gen[g] = gen.get(g, 0) + 1
        tot = sum(gen.values())
        top = sorted(gen.items(), key=lambda kv: -kv[1])[:6]
        print("      clade of %3d refs : %s" % (tot, ", ".join("%s x%d" % (g, n) for g, n in top)))
        shown += 1
        if shown >= 3: break
    if shown == 0: print("      (no reference leaves found in any enclosing clade)")
print("\n  Read the SMALLEST clade first. A query whose smallest reference-containing")
print("  clade is one genus is placed in that genus. A query whose smallest clade is")
print("  already mixed across families is unplaceable at genus level -- say so.")
PY
    else
      no "FastTree produced nothing -- see fasttree.err"
    fi
  else
    no "mafft produced nothing -- see mafft.err"
  fi
else
  no "skipped"
fi

# ================================================================ E
say "E. HOW TO DECIDE"
cat <<'EOT'
  THE CONTROLS COME FIRST, AND THEY CAN FAIL THIS TIME.

    CTRL1  LIQ_das_metabat_001   Cytobacillaceae  -> expect Bacilli, NOT a
                                 sulfate reducer
    CTRL2  LIQ_das_metabat_022   Desulfohalotomaculaceae -> the sharp control.
                                 If its reads ALSO come back "Peptococcaceae"
                                 or "Desulfotomaculum", then that is Kaiju's
                                 default answer for this whole family and the
                                 dsr contigs are placed no better than family.
    CTRL3  LIQ_das_vamb_017      Methanofollis -> must come back archaeal. If
                                 an archaeal contig returns sulfate reducers,
                                 the read-to-contig assignment is wrong and all
                                 of section A is void.

  THEN sections C and D, which do not depend on Kaiju at all. They are the
  standard method and they are the ones to quote. Agreement between an
  independent read-level signal and a protein phylogeny is a much stronger
  statement than either alone -- and disagreement is worth knowing before the
  chapter is written, not after.

  WATCH FOR THE REVERSE-TYPE dsr. Allochromatium, Thiobacillus and
  Desulfurivibrio references are in the set on purpose: rDsr belongs to sulfur
  OXIDISERS. If a query lands among them, the contig is not evidence of sulfate
  reduction and the chapter's central claim changes.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
