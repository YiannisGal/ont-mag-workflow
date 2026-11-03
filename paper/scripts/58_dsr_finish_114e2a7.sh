#!/usr/bin/env bash
# 58_dsr_finish.sh -- the two tests that close the dsrAB question.
#
#   bash scripts/run.sh scripts/58_dsr_finish_XXXXXXX.sh
#   SKIPMAP=1 bash scripts/run.sh scripts/58_dsr_finish_XXXXXXX.sh   (section B only)
#
# WHERE 57 LEFT IT. Three independent lines now point at LIQ_contig_19299
# belonging to ETO_LIQ_MAG_022 (LIQ_das_metabat_022, JBAIYS01):
#
#   tetranucleotide   nearest of all 29 genomes at 0.0157, inside that genome's
#                     own contig spread (p90 = 0.0161), and the runner-up is
#                     0.0285 -- 1.8x further away
#   read taxonomy     the same four taxa in the same rank order as MAG_022's
#                     own contig 2833 (Desulforadius, Desulfotomaculum,
#                     Desulfofalx, Bacillota bacterium)
#   no redundancy     none of the 26 contig proteins matches a MAG_022 protein
#                     at high identity over full length. Best dsr-related match
#                     is 22-23% -- superfamily level, not the same gene
#
# None of that is proof, and two things are still missing.
#
# COVERAGE (section A). If the contig belongs to that genome it was sequenced
# at that genome's depth. Step 56's per-contig read counts hinted at it
# (contig 2833 11.0 reads/kb, contig 19299 10.1) but those came from a
# 7-contig reference where reads had nowhere else to go. This maps the reads
# against ALL 29 genomes at once, so every read lands where it fits best and
# the depths are comparable.
#
# PHYLOGENY (section B). Volume XII asked for exactly this and step 56 could not
# do it because NCBI was unreachable -- but 57 section C found the real cause:
# wget fails on this machine while curl works. Same fetch, different client.
#
# Read-only except 58_dsr_finish/.

set -uo pipefail
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=..."; echo "STEP_STATUS: FAILED"; exit 1; }

OUT="$WS2/58_dsr_finish"; mkdir -p "$OUT"
export TMPDIR="${TMPDIR:-$WS2/tmp}"; mkdir -p "$TMPDIR"
exec > >(tee "$OUT/58.log") 2>&1
hr(){ printf '================================================================\n'; }
say(){ printf '\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---- %s\n' "$*"; }
ok(){ printf '  [OK]   %s\n' "$*"; }; no(){ printf '  [--]   %s\n' "$*"; }
bad(){ printf '  [FAIL] %s\n' "$*"; }
have(){ command -v "$1" >/dev/null 2>&1; }

echo "58 -- dsrAB, the closing tests"
echo "date : $(date -Is)"
echo "out  : $OUT"

PREV="$WS2/57_dsr_bin"
ASM="$WS2/04_polish/LIQ_asm_final.fasta"
READS="$WS2/00_input/liquid_reads.fastq.gz"
DSRFAA="$WS2/49_dsr/dsr_contig_proteins.faa"
CAND=LIQ_das_metabat_022
DSR_CONTIGS="LIQ_contig_12880 LIQ_contig_19299"
[ -s "$PREV/bin_paths.tsv" ] || { bad "run 57 first -- need $PREV/bin_paths.tsv"; echo "STEP_STATUS: FAILED"; exit 1; }

# ================================================================ A
say "A. COVERAGE, MEASURED COMPETITIVELY ACROSS ALL 29 GENOMES"

if [ "${SKIPMAP:-0}" = 1 ]; then
  no "SKIPMAP=1 -- section A skipped"
else
  REF="$OUT/all_genomes.fasta"
  if [ ! -s "$REF" ]; then
    : > "$REF"
    while IFS=$'\t' read -r B F T; do
      awk -v b="$B" '/^>/{ id=substr($1,2); sub(/.*~/,"",id); print ">" b "~" id; next } {print}' "$F" >> "$REF"
    done < "$PREV/bin_paths.tsv"
    samtools faidx "$ASM" $DSR_CONTIGS 2>/dev/null \
      | awk '/^>/{ id=substr($1,2); print ">UNBINNED~" id; next } {print}' >> "$REF"
  fi
  NSEQ=$(grep -c '^>' "$REF")
  SZ=$(awk '!/^>/{n+=length($0)} END{printf "%.1f", n/1e6}' "$REF")
  ok "reference: $NSEQ contigs, $SZ Mb (29 genomes + the 2 unbinned dsr contigs)"

  if [ ! -s "$OUT/all.bam" ]; then
    echo "  mapping $(basename "$READS") competitively -- this is the slow step"
    minimap2 -x map-ont -a --secondary=no -t 8 "$REF" "$READS" 2>"$OUT/minimap2.err" \
      | samtools sort -@ 4 -m 2G -o "$OUT/all.bam" - 2>>"$OUT/minimap2.err"
    samtools index "$OUT/all.bam" 2>/dev/null
  fi
  [ -s "$OUT/all.bam" ] || { bad "mapping failed -- see minimap2.err"; }

  if [ -s "$OUT/all.bam" ]; then
    samtools coverage "$OUT/all.bam" > "$OUT/coverage.tsv" 2>/dev/null
    ok "per-contig coverage -> ${OUT#$WS2/}/coverage.tsv"

    python3 - "$OUT/coverage.tsv" "$PREV/bin_paths.tsv" "$CAND" <<'PY'
import sys, statistics as st
cov, binsf, cand = sys.argv[1], sys.argv[2], sys.argv[3]
tax = {}
for line in open(binsf):
    b, p, t = line.rstrip('\n').split('\t')
    parts = [x for x in t.split(';') if len(x.split('__')) > 1 and x.split('__')[1]]
    tax[b] = ' / '.join(p2.split('__')[1] for p2 in parts[-2:]) if parts else t

per = {}
q = {}
for i, line in enumerate(open(cov)):
    f = line.rstrip('\n').split('\t')
    if i == 0 or f[0].startswith('#'): continue
    name = f[0]
    try:
        length = int(f[2]); depth = float(f[6]); cvg = float(f[5])
    except (ValueError, IndexError):
        continue
    b, _, c = name.partition('~')
    if b == 'UNBINNED':
        q[c] = (length, depth, cvg)
    elif length >= 5000:
        per.setdefault(b, []).append((c, length, depth, cvg))

def wmed(rows):
    ds = sorted(d for _, _, d, _ in rows)
    return st.median(ds) if ds else 0.0
def pct(rows, p):
    ds = sorted(d for _, _, d, _ in rows)
    if not ds: return 0.0
    i = min(len(ds)-1, max(0, int(round(p*(len(ds)-1)))))
    return ds[i]

print("\n  the two unbinned dsr contigs")
print("    %-20s %9s %10s %9s" % ("contig", "length", "meandepth", "covered%"))
for c, (L, d, cv) in q.items():
    print("    %-20s %9s %10.2f %8.1f%%" % (c, format(L, ','), d, cv))

print("\n  genome depth profiles (contigs >= 5 kb), ranked by closeness to contig 19299")
target = q.get('LIQ_contig_19299', (0, 0, 0))[1]
rows = []
for b, r in per.items():
    m = wmed(r)
    score = abs((m - target) / m) if m else 9e9
    rows.append((score, b, m, pct(r, 0.25), pct(r, 0.75), len(r)))
rows.sort()
print("    %-4s %-24s %9s %9s %9s %5s  %s" %
      ("rank", "genome", "med depth", "p25", "p75", "n", "taxonomy"))
for i, (s, b, m, p25, p75, n) in enumerate(rows[:8], 1):
    flag = "  <- contig 19299 sits inside this IQR" if p25 <= target <= p75 else ""
    print("    %-4d %-24s %9.2f %9.2f %9.2f %5d  %s%s" % (i, b, m, p25, p75, n, tax.get(b, ''), flag))
hit = [r for r in rows if r[1] == cand]
if hit and rows[0][1] != cand:
    s, b, m, p25, p75, n = hit[0]
    r = 1 + sum(1 for x in rows if x[0] < s)
    print("    %-4s %-24s %9.2f %9.2f %9.2f %5d  %s" % ("->", b + " (rank %d)" % r, m, p25, p75, n, tax.get(b, '')))

t2 = q.get('LIQ_contig_12880', (0, 0, 0))[1]
print("\n  and for contig 12880 (depth %.2f):" % t2)
rows2 = sorted(((abs((wmed(r) - t2) / wmed(r)) if wmed(r) else 9e9, b, wmed(r)) for b, r in per.items()))
for i, (s, b, m) in enumerate(rows2[:4], 1):
    print("    %-4d %-24s med depth %8.2f   %s" % (i, b, m, tax.get(b, '')))

print("""
  WHAT THIS CAN AND CANNOT SAY. Matching depth is consistent with belonging to
  a genome; it is not diagnostic, because two organisms at similar abundance
  have similar depth. Depth that does NOT match is the informative direction --
  it would rule the candidate out. Read this as a test the hypothesis had to
  survive, not as the thing that proves it.""")
PY
  fi
fi

# ================================================================ B
say "B. DsrAB PHYLOGENY -- the method Volume XII actually asked for"

EUT="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CURL="curl -s --max-time 60"
fetch(){  # $1 outfile  $2 retmax  $3 term
  local out="$1" n="$2" term="$3" ids
  ids=$($CURL --get --data-urlencode "term=$term" \
        "$EUT/esearch.fcgi?db=protein&retmax=$n" 2>/dev/null \
        | grep -oP '(?<=<Id>)[0-9]+(?=</Id>)' | tr '\n' ',' | sed 's/,$//')
  if [ -z "$ids" ]; then printf '      %-56s %s\n' "$term" "0"; return 1; fi
  $CURL "$EUT/efetch.fcgi?db=protein&id=${ids}&rettype=fasta&retmode=text" >> "$out" 2>/dev/null
  printf '      %-56s %s\n' "$term" "$(printf '%s' "$ids" | tr ',' '\n' | wc -l)"
  sleep 0.4
}

if ! $CURL -o /dev/null "$EUT/einfo.fcgi"; then
  bad "NCBI still unreachable with curl. Section B cannot run."
else
  ok "NCBI reachable via curl"
  if [ ! -s "$OUT/refs_raw.faa" ]; then
    : > "$OUT/refs_raw.faa"
    for SUB in alpha beta; do
      sub "dissimilatory sulfite reductase $SUB subunit"
      fetch "$OUT/refs_raw.faa" 90 "\"dissimilatory sulfite reductase $SUB subunit\"[Protein Name] AND refseq[filter]"
      # reductive lineages, then the reverse-type sulfur OXIDISERS on purpose
      for ORG in Desulfosporosinus Desulfotomaculum Desulforamulus Desulfofalx \
                 Pelotomaculum Desulfitobacterium Desulfovibrio Desulfobacterium \
                 Desulfobacter Archaeoglobus Thermodesulfovibrio \
                 Desulfurivibrio Allochromatium Thiobacillus Chlorobium; do
        fetch "$OUT/refs_raw.faa" 5 "\"dissimilatory sulfite reductase $SUB subunit\"[Protein Name] AND \"$ORG\"[Organism]"
      done
    done
  fi

  if [ -s "$OUT/refs_raw.faa" ]; then
    awk '
      /^>/ { acc=substr($1,2); org=""
             if (match($0,/\[[^]]*\]$/)) org=substr($0,RSTART+1,RLENGTH-2)
             gsub(/[^A-Za-z0-9]+/,"_",org); if(org=="") org="unknown"
             if (acc in seen) { skip=1; next }
             seen[acc]=1; skip=0; print ">REF|" org "|" acc; next }
      { if(!skip) print }' "$OUT/refs_raw.faa" > "$OUT/refs.faa"
    NR_=$(grep -c '^>' "$OUT/refs.faa" | head -1)
    ok "${NR_:-0} unique reference proteins"
    grep '^>' "$OUT/refs.faa" | awk -F'|' '{split($2,a,"_"); print a[1]}' \
      | sort | uniq -c | sort -rn | head -18 | sed 's/^/    /'
  fi
fi

QIDS="LIQ_contig_12880_3 LIQ_contig_12880_2 LIQ_contig_19299_15 LIQ_contig_19299_14 LIQ_contig_19299_13"
awk -v ids="$QIDS" '
  BEGIN{ n=split(ids,a," "); for(i=1;i<=n;i++) want[a[i]]=1 }
  /^>/{ id=substr($1,2); p=(id in want); if(p) print ">QUERY|" id; next }
  { if(p) print }' "$DSRFAA" > "$OUT/queries.faa"
NQ=$(grep -c '^>' "$OUT/queries.faa" | head -1); NQ=${NQ:-0}
ok "$NQ query proteins"

# ================================================================ C
say "C. NEAREST NAMED RELATIVES"
if [ -s "$OUT/refs.faa" ] && [ "$NQ" -gt 0 ]; then
  makeblastdb -in "$OUT/refs.faa" -dbtype prot -out "$OUT/refdb" >/dev/null 2>&1
  blastp -query "$OUT/queries.faa" -db "$OUT/refdb" \
         -outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' \
         -max_target_seqs 6 -max_hsps 1 -num_threads 8 > "$OUT/blastp.tsv" 2>/dev/null
  awk -F'\t' '
    { if($1!=q){ q=$1; printf "\n  %s\n", q
                 printf "      %-46s %6s %7s %6s %6s\n","nearest named relative","pid%","aln","qlen","slen" }
      split($2,s,"|"); printf "      %-46s %6.1f %7d %6d %6d\n", s[2]" ("s[3]")", $3, $4, $5, $6 }' "$OUT/blastp.tsv"
else
  no "skipped"
fi

# ================================================================ D
say "D. PLACEMENT"
if [ -s "$OUT/refs.faa" ] && [ "$NQ" -gt 0 ]; then
  cat "$OUT/queries.faa" "$OUT/refs.faa" > "$OUT/all.faa"
  mafft --auto --quiet --thread 8 "$OUT/all.faa" > "$OUT/all.aln" 2>"$OUT/mafft.err"
  if [ -s "$OUT/all.aln" ]; then
    FastTree -lg -quiet "$OUT/all.aln" > "$OUT/dsr.tree" 2>"$OUT/fasttree.err"
  fi
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
OXID = {'Allochromatium', 'Thiobacillus', 'Chlorobium', 'Desulfurivibrio'}
for qy in [l for l in leaves(root) if l.startswith('QUERY')]:
    print("\n  %s" % qy)
    path = []
    def go(n):
        if not n['ch']: return n['lab'] == qy
        if any(go(c) for c in n['ch']): path.append(n); return True
        return False
    go(root)
    shown = 0
    for n in path:
        refs = [x for x in leaves(n) if x.startswith('REF') and x != qy]
        if not refs: continue
        gen = {}
        for r in refs:
            g = r.split('|')[1].split('_')[0]
            gen[g] = gen.get(g, 0) + 1
        top = sorted(gen.items(), key=lambda kv: -kv[1])[:6]
        ox = sum(v for g, v in gen.items() if g in OXID)
        note = "   <-- REVERSE-TYPE dsr present" if ox and ox >= 0.3 * sum(gen.values()) else ""
        print("      clade of %3d refs : %s%s" %
              (sum(gen.values()), ", ".join("%s x%d" % (g, n2) for g, n2 in top), note))
        shown += 1
        if shown >= 3: break
    if shown == 0: print("      (no reference leaves in any enclosing clade)")
print("""
  Read the SMALLEST clade first. If it is dominated by Desulfotomaculum /
  Desulforamulus / Desulfofalx, the contig is a Desulfotomaculales sulfate
  reducer and consistent with MAG_022. If it is Desulfosporosinus, that is the
  12880 organism. If a query sits with Allochromatium / Thiobacillus /
  Chlorobium, it is reverse-type dsr from a sulfur OXIDISER and does not
  support sulfate reduction at all.""")
PY
  else
    no "no tree -- see mafft.err / fasttree.err"
  fi
else
  no "skipped"
fi

say "E. THE DECISION THIS RUN SUPPORTS"
cat <<'EOT'
  For LIQ_contig_19299, four things now have to agree before it is written as
  MAG_022's operon: composition (done, compatible), read taxonomy (done,
  matches), no protein redundancy (done, none), depth (section A) and phylogeny
  (section D). If depth is inside MAG_022's contig IQR and the DsrAB places in
  Desulfotomaculales, that is as far as short-read-free evidence goes without
  re-binning, and the chapter should say:

      "The dsrAB operon on LIQ_contig_19299 is compositionally, taxonomically
       and in coverage consistent with ETO_LIQ_MAG_022 (JBAIYS01,
       Desulfohalotomaculaceae), the genome that carries sat, aprAB and dsrMK
       but no dsrAB. It is most parsimoniously that genome's missing operon,
       left unbinned at 14x coverage on a 20.7 kb contig."

  with "most parsimoniously" doing real work -- it is an attribution, not an
  identification.

  If depth or phylogeny disagrees, MAG_022 keeps its gap and the lagoon has a
  second, distinct Desulfotomaculales sulfate reducer that never binned. Say
  that instead. Both are publishable; only one is true.

  LIQ_contig_12880 is separate either way: Desulfosporosinus by reads, nearest
  no MAG in composition, and its own organism.
EOT

echo
echo "STEP_STATUS: COMPLETE $(date -Is)"
