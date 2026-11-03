#!/usr/bin/env bash
# 57_dsr_binning_test.sh -- does LIQ_contig_19299 belong to ETO_LIQ_MAG_022?
#
#   bash scripts/run.sh scripts/57_dsr_binning_test_XXXXXXX.sh
#
# WHAT 56 ESTABLISHED. The controls passed and calibrated the method:
#   CTRL1 Cytobacillaceae contig -> Bacillaceae/Robertmurraya/Lysinibacillus
#   CTRL3 Methanofollis contigs  -> Methanofollis formosanus, 9,264 reads
#   CTRL2 Desulfohalotomaculaceae (LIQ_das_metabat_022, JBAIYS01), contig 2833
#         -> Desulforadius tongensis 121, Desulfofalx alkaliphila 94,
#            Desulfotomaculum sp. 92, Bacillota bacterium 72
#
# And LIQ_contig_19299 -> Bacillota bacterium 34, Desulforadius tongensis 32,
# Desulfotomaculum sp. 19, Desulfofalx alkaliphila 18.
#
# THE SAME FOUR TAXA IN THE SAME ORDER as the known JBAIYS01 genome. Two
# readings, and they are not equally interesting:
#   (i)  Kaiju cannot resolve inside Desulfotomaculales, so this is its generic
#        answer and contig 19299 is merely "some Desulfotomaculales".
#   (ii) Contig 19299 IS FROM MAG_022 -- the genome that carries sat, aprAB and
#        dsrMK and was recorded as lacking dsrAB. If so, the chapter gains a
#        named, near-complete sulfate reducer and loses its central absence.
#
# Volume XII withdrew the MAG_022 proposal on GC alone. GC is one number; this
# script uses the signal binners actually use.
#
# LIQ_contig_12880 is a separate organism and the same tests apply to it. Its
# Desulfosporosinus signature appears in NEITHER control, so it is not a
# generic answer.
#
# Read-only except 57_dsr_bin/.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/57_dsr_bin"; mkdir -p "$OUT"
export TMPDIR="${TMPDIR:-$WS2/tmp}"; mkdir -p "$TMPDIR"
exec > >(tee "$OUT/57.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }; no(){ printf '  [--]   %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
have(){ command -v "$1" >/dev/null 2>&1; }
nlines(){ grep -c '^>' "$1" 2>/dev/null || true; }

echo "57 -- dsr contig binning test"
echo "date : $(date -Is)"
echo "out  : $OUT"

ASM="$WS2/04_polish/LIQ_asm_final.fasta"
DSRFAA="$WS2/49_dsr/dsr_contig_proteins.faa"
BINS="$WS2/06_binqc/union_bins"
PROT="$WS2/06_binqc/checkm2/protein_files"
GTDB_B="$WS2/08_taxonomy/gtdbtk/classify/gtdbtk.bac120.summary.tsv"
GTDB_A="$WS2/08_taxonomy/gtdbtk/classify/gtdbtk.ar53.summary.tsv"
DSR_CONTIGS="LIQ_contig_12880 LIQ_contig_19299"
CAND=LIQ_das_metabat_022

# ================================================================ A
say "A. THE FINAL GENOME SET, AND THE QUERY CONTIGS"

awk -F'\t' 'FNR>1 && $1!=""{print $1"\t"$2}' "$GTDB_B" "$GTDB_A" 2>/dev/null | sort -u > "$OUT/mag_tax.tsv"
NM=$(wc -l < "$OUT/mag_tax.tsv")
ok "$NM genomes classified by GTDB-Tk"
[ "$NM" -lt 20 ] && { bad "expected ~29 -- check the summary paths"; }

: > "$OUT/bin_paths.tsv"
while IFS=$'\t' read -r B T; do
  F="$(find "$BINS" -maxdepth 1 -name "${B}.fa" -o -maxdepth 1 -name "${B}.fasta" 2>/dev/null | head -1)"
  [ -n "$F" ] && printf '%s\t%s\t%s\n' "$B" "$F" "$T" >> "$OUT/bin_paths.tsv"
done < "$OUT/mag_tax.tsv"
ok "$(wc -l < "$OUT/bin_paths.tsv") bin fastas located"

samtools faidx "$ASM" $DSR_CONTIGS > "$OUT/dsr_contigs.fasta" 2>/dev/null
ok "query contigs: $(nlines "$OUT/dsr_contigs.fasta")"

sub "the candidate, for reference"
awk -F'\t' -v b="$CAND" '$1==b{print "  "$1"\n  "$3}' "$OUT/bin_paths.tsv"

# ================================================================ B
say "B. QUERY PROTEINS -- extracted properly this time"
# 56 failed here: seqkit grep -n needs the FULL header including the prodigal
# description, and the awk fallback never ran because a broken grep -c gave
# NQ the two-line value "0\n0". Plain awk, counted with wc.
QIDS="LIQ_contig_12880_3 LIQ_contig_12880_2 LIQ_contig_19299_15 LIQ_contig_19299_14 LIQ_contig_19299_13"
awk -v ids="$QIDS" '
  BEGIN{ n=split(ids,a," "); for(i=1;i<=n;i++) want[a[i]]=1 }
  /^>/{ id=substr($1,2); p=(id in want); if(p) print ">QUERY|" id; next }
  { if(p) print }' "$DSRFAA" > "$OUT/queries.faa"
NQ=$(grep -c '^>' "$OUT/queries.faa" | head -1)
NQ=${NQ:-0}
ok "$NQ query proteins extracted"
awk '/^>/{if(n)printf " %d aa\n", n; printf "    %-32s", substr($0,2); n=0; next}{n+=length($0)}
     END{if(n)printf " %d aa\n", n}' "$OUT/queries.faa"
echo "    reference lengths: DsrA ~430 aa, DsrB ~360 aa"

# ================================================================ C
say "C. IS NCBI REACHABLE AT ALL? (56 sections B-D depend on it)"
for probe in "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi" \
             "https://rest.uniprot.org/" "https://www.google.com"; do
  if wget -q --spider --timeout=10 --tries=1 "$probe" 2>/dev/null; then
    printf '  [OK]   wget  %s\n' "$probe"
  elif have curl && curl -s -o /dev/null -m 10 "$probe" 2>/dev/null; then
    printf '  [OK]   curl  %s   (wget failed, curl works)\n' "$probe"
  else
    printf '  [--]   unreachable: %s\n' "$probe"
  fi
done
echo "  proxy env: http_proxy=${http_proxy:-unset} https_proxy=${https_proxy:-unset}"
echo "  (if curl works where wget does not, 56 can be re-run with curl)"

# ================================================================ D
say "D. TETRANUCLEOTIDE COMPOSITION -- the signal binners actually use"
python3 - "$OUT" "$OUT/bin_paths.tsv" "$OUT/dsr_contigs.fasta" <<'PY'
import sys, os, math, itertools

OUT, BINS_TSV, QFA = sys.argv[1], sys.argv[2], sys.argv[3]

def rc(s):
    c = {'A':'T','C':'G','G':'C','T':'A'}
    return ''.join(c.get(b,'N') for b in reversed(s))

KM, seen = [], set()
for t in itertools.product('ACGT', repeat=4):
    k = ''.join(t)
    if k in seen: continue
    seen.add(k); seen.add(rc(k)); KM.append(k)
IDX = {}
for i,k in enumerate(KM):
    IDX[k] = i; IDX[rc(k)] = i

def read_fasta(p):
    name, buf = None, []
    with open(p) as fh:
        for line in fh:
            if line.startswith('>'):
                if name: yield name, ''.join(buf)
                name = line[1:].split()[0]; buf = []
            else: buf.append(line.strip())
    if name: yield name, ''.join(buf)

def prof(seq):
    v = [0.0]*len(KM); n = 0
    seq = seq.upper()
    for i in range(len(seq)-3):
        j = IDX.get(seq[i:i+4])
        if j is not None:
            v[j] += 1; n += 1
    return [x/n for x in v] if n else None

def gc(seq):
    seq = seq.upper(); g = seq.count('G')+seq.count('C')
    at = g + seq.count('A')+seq.count('T')
    return 100.0*g/at if at else 0.0

def dist(a,b):
    return math.sqrt(sum((x-y)**2 for x,y in zip(a,b)))

queries = []
for nm, sq in read_fasta(QFA):
    queries.append((nm, len(sq), gc(sq), prof(sq)))

bins = []
for line in open(BINS_TSV):
    b, path, tax = line.rstrip('\n').split('\t')
    seqs = list(read_fasta(path))
    whole = ''.join(s for _, s in seqs)
    if len(whole) < 50000: continue
    wp = prof(whole)
    # within-genome spread: contigs >= 5 kb against their own genome profile
    spread = sorted(dist(prof(s), wp) for _, s in seqs if len(s) >= 5000)
    bins.append(dict(name=b, tax=tax, prof=wp, gc=gc(whole), n=len(seqs),
                     size=len(whole), spread=spread))

def q(a, p):
    if not a: return float('nan')
    i = min(len(a)-1, max(0, int(round(p*(len(a)-1)))))
    return a[i]

def short(t):
    parts = [x for x in t.split(';') if len(x.split('__')) > 1 and x.split('__')[1]]
    return ' / '.join(p.split('__')[1] for p in parts[-2:]) if parts else t

for nm, L, g, p in queries:
    print("\n  %s   %s bp   GC %.2f%%" % (nm, format(L, ','), g))
    rows = sorted(((dist(p, b['prof']), b) for b in bins), key=lambda r: r[0])
    print("    %-4s %-24s %9s %8s %10s %10s  %s" %
          ("rank", "MAG", "TNF dist", "dGC", "own p50", "own p90", "taxonomy"))
    for i, (d, b) in enumerate(rows[:8], 1):
        flag = ""
        if b['spread']:
            if d <= q(b['spread'], 0.90): flag = "  <= its own p90  COMPATIBLE"
            elif d <= max(b['spread']):   flag = "  <= its own max"
        print("    %-4d %-24s %9.4f %+8.2f %10s %10s  %s%s" %
              (i, b['name'], d, g-b['gc'],
               ("%.4f" % q(b['spread'],0.50)) if b['spread'] else "-",
               ("%.4f" % q(b['spread'],0.90)) if b['spread'] else "-",
               short(b['tax']), flag))
    cand = [b for b in bins if b['name'] == 'LIQ_das_metabat_022']
    if cand and rows[0][1]['name'] != 'LIQ_das_metabat_022':
        b = cand[0]; d = dist(p, b['prof'])
        r = 1 + sum(1 for dd,_ in rows if dd < d)
        print("    %-4s %-24s %9.4f %+8.2f %10s %10s  %s" %
              ("->", b['name']+" (rank %d)"%r, d, g-b['gc'],
               ("%.4f" % q(b['spread'],0.50)) if b['spread'] else "-",
               ("%.4f" % q(b['spread'],0.90)) if b['spread'] else "-",
               short(b['tax'])))

print("""
  HOW TO READ THIS. "own p50/p90" is how far that genome's OWN contigs sit from
  its own whole-genome profile. It is the yardstick: a query closer to a genome
  than that genome's own contigs typically are is compositionally compatible
  with belonging to it. A query further away than the genome's worst contig is
  not. Compatibility is necessary, not sufficient -- related genomes have
  similar composition, which is exactly why section E exists.""")
PY

# ================================================================ E
say "E. SEQUENCE-LEVEL: is either contig already present in a MAG?"
cat $(awk -F'\t' '{print $2}' "$OUT/bin_paths.tsv") > "$OUT/all_mags.fna" 2>/dev/null
ok "concatenated MAG contigs: $(nlines "$OUT/all_mags.fna")"
if have blastn; then
  makeblastdb -in "$OUT/all_mags.fna" -dbtype nucl -out "$OUT/magdb" >/dev/null 2>&1
  blastn -query "$OUT/dsr_contigs.fasta" -db "$OUT/magdb" -outfmt '6 qseqid sseqid pident length qstart qend evalue' \
         -max_target_seqs 10 -num_threads 8 -evalue 1e-20 > "$OUT/blastn.tsv" 2>/dev/null
  if [ -s "$OUT/blastn.tsv" ]; then
    awk -F'\t' '$4>=500{printf "    %-18s -> %-26s %6.2f%% id  %6d bp  q%d-%d\n",$1,$2,$3,$4,$5,$6}' "$OUT/blastn.tsv" | head -20
    [ -z "$(awk -F'\t' '$4>=500' "$OUT/blastn.tsv")" ] && echo "    no alignment >= 500 bp to any binned contig"
  else
    echo "    no hits at all -- neither contig is a duplicate of binned sequence"
  fi
else
  no "blastn not available"
fi

sub "and are the contig proteins ALREADY in the candidate genome?"
CP="$(find "$PROT" -name "${CAND}*.faa" | head -1)"
if [ -n "$CP" ] && [ -s "$DSRFAA" ] && have blastp; then
  ok "candidate proteome: ${CP#$WS2/}  ($(nlines "$CP") proteins)"
  makeblastdb -in "$CP" -dbtype prot -out "$OUT/canddb" >/dev/null 2>&1
  blastp -query "$DSRFAA" -db "$OUT/canddb" -outfmt '6 qseqid sseqid pident length qlen slen evalue' \
         -max_target_seqs 1 -num_threads 8 > "$OUT/blastp_cand.tsv" 2>/dev/null
  NH=$(wc -l < "$OUT/blastp_cand.tsv")
  echo "    $NH of 26 contig proteins have a hit in $CAND"
  awk -F'\t' '{printf "    %-24s -> %-22s %6.2f%% id  aln %4d / q%4d\n",$1,$2,$3,$4,$5}' "$OUT/blastp_cand.tsv" | head -26
  echo
  echo "    A contig that BELONGS to this genome should show few or no"
  echo "    high-identity full-length matches -- a genome does not carry the"
  echo "    same gene twice. Many >90% full-length matches mean the contig is a"
  echo "    REDUNDANT copy, i.e. a strain variant or a different organism."
else
  no "candidate proteome or blastp unavailable"
fi

# ================================================================ F
say "F. INDEPENDENT TAXONOMY -- kraken2 on the contigs themselves"
KDB="$HOME/databases/kraken2_std16"
if have kraken2 && [ -s "$KDB/hash.k2d" ]; then
  ok "kraken2 + $KDB"
  CTRLC="LIQ_contig_2833 LIQ_contig_1679"
  samtools faidx "$ASM" $DSR_CONTIGS $CTRLC > "$OUT/k2_input.fasta" 2>/dev/null
  kraken2 --db "$KDB" --threads 8 --report "$OUT/k2.report" --use-names \
          "$OUT/k2_input.fasta" > "$OUT/k2.out" 2>"$OUT/k2.err"
  if [ -s "$OUT/k2.out" ]; then
    awk -F'\t' '{printf "    %-4s %-20s %10s  %s\n",$1,$2,$4,$3}' "$OUT/k2.out"
    echo
    echo "    LIQ_contig_2833 is a MAG_022 contig and LIQ_contig_1679 a"
    echo "    Methanofollis contig -- both are here as controls. kraken2 uses a"
    echo "    different database and algorithm from kaiju, so agreement between"
    echo "    them is real corroboration rather than the same result twice."
  else
    no "kraken2 produced no output -- see k2.err"
  fi
else
  no "kraken2 not on PATH, or $KDB/hash.k2d missing"
  echo "         (names.dmp is there, so the db directory may be taxonomy-only)"
fi

say "G. WHAT EACH OUTCOME MEANS"
cat <<'EOT'
  FOR LIQ_contig_19299:

  D says COMPATIBLE with LIQ_das_metabat_022, E finds no redundant proteins,
  F calls it Desulfotomaculales
      -> report that the dsrAB operon most likely belongs to ETO_LIQ_MAG_022
         and failed to bin, and say plainly that it is a composition-plus-
         taxonomy argument, not proof. The chapter then has a named sulfate
         reducer, and "no MAG carries dsrAB" becomes "MAG_022's dsrAB did not
         bin".

  D ranks it far from metabat_022, or E finds many high-identity full-length
  protein matches
      -> it is a DIFFERENT Desulfotomaculales organism. MAG_022 keeps its
         missing dsrAB and the chapter keeps the absence, now bounded: a second
         sulfate reducer is present in the assembly and was never binned.

  FOR LIQ_contig_12880: Desulfosporosinus appears in no control, so the read
  signal is already specific. If D puts it far from every MAG, it is simply an
  unbinned organism with no counterpart in the catalogue -- which is a clean,
  reportable statement and matches its 33 reads of coverage.

  A NOTE ON WHAT COVERAGE ALLOWS. 33 reads on 6.4 kb and 209 on 20.7 kb is
  roughly 7x and 14x. No binner places contigs that short at that depth with
  confidence, so "not binned" needs no further explanation than this.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
