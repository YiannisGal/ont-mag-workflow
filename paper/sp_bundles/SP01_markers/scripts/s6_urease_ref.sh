#!/usr/bin/env bash
# SP01 step 5a -- ureC (K01428) read-level reference set, built the step-34 way:
#   1. verified assembly proteins: SED/LIQ assembly proteins (34_unbinned/prot) passing the K01428 curated threshold
#   2. UniProt gene:ureC, kept only if the protein name contains "urease" (step-34 name-collision filter)
#   3. DECOY amidohydrolase relatives (reviewed): a read whose best hit is a decoy is NOT counted as ureC
#   4. CONTROLS rpoB / glnA / tuf: the identical CTRL sequences used by step 34 (copied from its marker_ref.faa)
set -euo pipefail
U=$OUT/urease; mkdir -p "$U"; cd "$U"
S34=$WS2/34_unbinned
cp "$PROFILES/K01428.hmm" K01428.hmm
for s in SED LIQ; do
  $HMMSEARCH --cpu "$CPU" --tblout "asm_${s}.tbl" K01428.hmm "$S34/prot/${s}_all.faa" > /dev/null
done
thr=$(awk -F'\t' '$1=="K01428"{print $2" "$3}' "$KOLIST")
echo "K01428 threshold/score_type: $thr"
: > urease_ref.faa
$PY - "$thr" "$S34" <<'PY' >> urease_ref.faa
import sys
t, st = sys.argv[1].split(); t = float(t); S34 = sys.argv[2]
for s in ('SED', 'LIQ'):
    want = set()
    for l in open(f'asm_{s}.tbl'):
        if l.startswith('#'): continue
        f = l.split(); sc = float(f[8]) if st == 'domain' else float(f[5])
        if sc >= t: want.add(f[0])
    sys.stderr.write(f'  verified assembly ureC {s}: {len(want)}\n')
    keep = False
    for l in open(f'{S34}/prot/{s}_all.faa'):
        if l.startswith('>'):
            n = l[1:].split()[0]; keep = n in want
            if keep: print(f'>ureC|ASM_{s}_{n}')
        elif keep: print(l.rstrip())
PY
fetch(){  # fetch <label> <query> <keyword or ->
  local lab=$1 q=$2 kw=$3
  curl -sf -m 120 "https://rest.uniprot.org/uniprotkb/search?format=fasta&size=500&query=$q" -o ".up_$lab.faa"
  $PY - "$lab" "$kw" ".up_$lab.faa" <<'PY' >> urease_ref.faa
import sys
lab, kw, p = sys.argv[1:4]; keep = False; k = d = 0
for l in open(p):
    if l.startswith('>'):
        keep = kw == '-' or kw in l.lower()
        k += keep; d += not keep
        if keep: print(f'>{lab}|{l[1:].rstrip()}')
    elif keep: print(l.rstrip())
sys.stderr.write(f'  {lab}: kept {k}, name-filter dropped {d}\n')
PY
  rm -f ".up_$lab.faa"
}
fetch ureC        '%28gene%3AureC%29'                                   urease
fetch DECOYpyrC   '%28gene%3ApyrC%29+AND+%28reviewed%3Atrue%29'         -
fetch DECOYhutI   '%28gene%3AhutI%29+AND+%28reviewed%3Atrue%29'         -
fetch DECOYallB   '%28gene%3AallB%29+AND+%28reviewed%3Atrue%29'         -
# controls: exactly step 34's
$PY - "$S34/marker_ref.faa" <<'PY' >> urease_ref.faa
import sys; keep = False
for l in open(sys.argv[1]):
    if l.startswith('>'): keep = l[1:].split('|')[0] in ('CTRLrpoB', 'CTRLglnA', 'CTRLtuf')
    if keep: print(l.rstrip())
PY
grep '^>' urease_ref.faa | cut -d'|' -f1 | sort | uniq -c
$DIAMOND makedb --in urease_ref.faa -d urease_ref --quiet
