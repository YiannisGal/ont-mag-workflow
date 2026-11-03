#!/usr/bin/env bash
# =============================================================================
# 15_marker_trees.sh  --  resolve the three annotation traps by phylogeny
#
# THREE CLAIMS CURRENTLY REST ON UNVERIFIED HMM HITS
#
#   1. mcrA  (7 MAGs)  "these are methanogens"
#      Alkyl-CoM reductases (ACR) of anaerobic alkane oxidisers — Ca.
#      Methanoliparia, Syntrophoarchaeum, some Archaeoglobi and Bathyarchaeia —
#      are homologous to MCR and score identically against K00399. workspace_1
#      had roughly 50% false positives on "complete mcrABG". Until MCR and ACR
#      are separated on a tree, a complete mcr set is a CANDIDATE.
#
#      Note the supporting evidence already in hand: all seven sit in bona fide
#      methanogen families by GTDB (Methanosarcinaceae, Methanofollaceae,
#      Methanomicrobiaceae, Methanobacteriaceae) and carry complete mtr operons
#      and H4MPT pathways, which ACR-only lineages do not. The tree is the
#      formal confirmation of what the pathway evidence already suggests.
#
#   2. nifH  (14 MAGs, incl. ALL 8 methanogens)  "nitrogen fixation"
#      PREDICTION TO TEST: most archaeal hits are cfbC/cfbD, the coenzyme F430
#      biosynthesis proteins, which are nitrogenase homologues. Every methanogen
#      needs F430 for MCR, so every methanogen carries them. bchL/chlL
#      (protochlorophyllide reductase) and nflH are the other usual suspects.
#      If the archaeal sequences branch with cfbC/cfbD, there is NO nitrogen
#      fixation claim to make — which is how that error reaches the literature.
#
#   3. mttB  (3 Methanosarcinaceae)  "methylamine methanogenesis"
#      Already corroborated independently: those three MAGs, and only those,
#      carry pyrrolysine biosynthesis (pyl) operons, and pyrrolysine exists
#      essentially only to build these enzymes. The tree is confirmation of a
#      conclusion that already has independent support, not the sole evidence.
#
# METHOD
#   Reference proteins are fetched from NCBI by DESCRIPTION rather than from a
#   hardcoded accession list, so the sets stay current and are reproducible from
#   the query strings recorded here. MAFFT-linsi -> trimAl -> IQ-TREE with
#   ModelFinder, 1000 UFBoot and SH-aLRT.
#
# INPUT    11_kofam/verify/*.faa   (extracted at step 11b)
# OUTPUT   15_marker_trees/<marker>/  tree + per-MAG verdict table
#
# RUNTIME  20-50 min.
# USAGE    ws2 && bash scripts/15_marker_trees.sh 2>&1 | tee logs/15_markers.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="15_marker_trees"
export TMPDIR="${TMPDIR}"; mkdir -p "$TMPDIR"
V="$WS2/11_kofam/verify"; M="$WS2/15_marker_trees"; mkdir -p "$M"

echo "######## 15 — MARKER GENE TREES ########"
echo "  started: $(date -Is)"
runlog "15_marker_trees started"
[ -d "$V" ] || { echo "  !! $V missing — run 11b first"; exit 1; }

EUT="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
fetch() {  # fetch <outfile> <retmax> <query> [fallback query ...]
  # Tries each query in turn and accumulates. NCBI phrasing is not guessable:
  # RefSeq calls it "methyl-coenzyme M reductase subunit alpha", NOT
  # "... alpha subunit", and the first version of this script fetched ZERO mcrA
  # references because of that. Always supply fallbacks and check the count.
  local out="$1" n="$2"; shift 2
  [ -s "$out" ] && { echo "      cached: $(grep -c '^>' "$out") seqs"; return 0; }
  : > "$out.tmp"
  local q enc ids
  for q in "$@"; do
    enc=$(printf '%s' "$q" | sed 's/ /+/g; s/"/%22/g; s/\[/%5B/g; s/\]/%5D/g; s/(/%28/g; s/)/%29/g')
    ids=$(wget -qO- "$EUT/esearch.fcgi?db=protein&term=${enc}&retmax=${n}" \
          | grep -oE '<Id>[0-9]+</Id>' | grep -oE '[0-9]+' | paste -sd,)
    if [ -z "$ids" ]; then echo "      (0)  $q"; sleep 1; continue; fi
    wget -qO- "$EUT/efetch.fcgi?db=protein&id=${ids}&rettype=fasta&retmode=text" >> "$out.tmp"
    echo "      (+$(grep -c '^>' "$out.tmp" 2>/dev/null || echo 0) cumulative)  $q"
    sleep 1
    [ "$(grep -c '^>' "$out.tmp" 2>/dev/null || echo 0)" -ge "$n" ] && break
  done
  if [ "$(grep -c '^>' "$out.tmp" 2>/dev/null || echo 0)" -eq 0 ]; then
    rm -f "$out.tmp"; echo "      !! NO references retrieved for this group"; return 1
  fi
  run_tool seqkit rmdup -s "$out.tmp" -o "$out" 2>/dev/null || mv "$out.tmp" "$out"
  rm -f "$out.tmp"
  echo "      fetched: $(grep -c '^>' "$out") unique seqs"
}

label_refs() {  # label_refs <fasta> <tag>  -> prefix every header with the group
  awk -v t="$2" '/^>/{sub(/^>/,""); split($0,a," "); print ">"t"|"a[1]; next}{print}' "$1"
}

build_marker() {  # build_marker <name> <candidates.faa>
  local mk="$1" cand="$2"
  [ -s "$cand" ] || { echo "  $mk: no candidate sequences — skipping"; return 0; }
  local d="$M/$mk"; mkdir -p "$d/refs"
  echo
  echo "=================================================================="
  echo "  $mk  —  $(grep -c '^>' "$cand") candidate sequences from our MAGs"
  echo "=================================================================="

  local all="$d/${mk}_all.faa"
  if [ ! -s "$all" ]; then
    : > "$all"
    # our sequences, tagged
    awk '/^>/{sub(/^>/,""); print ">QUERY|"$0; next}{print}' "$cand" >> "$all"
    case "$mk" in
      mcrA)
        echo "    fetching references ..."
        # true MCR from cultured methanogens
        fetch "$d/refs/mcr.faa" 60 \
          'mcrA[Gene Name] AND refseq[filter]' \
          '"methyl-coenzyme M reductase subunit alpha"[Title]' \
          '"methyl-coenzyme M reductase"[Title] AND alpha[Title]' \
          '"coenzyme-B sulfoethylthiotransferase"[Title]'
        # ACR / divergent MCR from the alkane-oxidising lineages — the exact
        # sequences that make an mcrA hit ambiguous
        fetch "$d/refs/acr.faa" 40 \
          '"alkyl-coenzyme M reductase"[Title]' \
          '"methyl-coenzyme M reductase"[Title] AND (Syntrophoarchaeum[Organism] OR Methanoliparia[Organism] OR "Candidatus Methanoliparum"[Organism] OR Archaeoglobi[Organism] OR Bathyarchaeia[Organism])' \
          'mcrA[Gene Name] AND (Syntrophoarchaeum[Organism] OR "Candidatus Methanoliparum"[Organism])'
        [ -s "$d/refs/mcr.faa" ] && label_refs "$d/refs/mcr.faa" MCR >> "$all"
        [ -s "$d/refs/acr.faa" ] && label_refs "$d/refs/acr.faa" ACR >> "$all"
        ;;
      nifH)
        echo "    fetching references ..."
        fetch "$d/refs/nifh.faa" 50 '"nitrogenase iron protein"[Protein Name] AND refseq[filter]'
        fetch "$d/refs/cfb.faa"  30 '"coenzyme F430 synthetase"[All Fields] OR CfbC[All Fields]'
        fetch "$d/refs/bchl.faa" 30 '"light-independent protochlorophyllide reductase iron-sulfur ATP-binding protein"[Protein Name]'
        [ -s "$d/refs/nifh.faa" ] && label_refs "$d/refs/nifh.faa" NIFH   >> "$all"
        [ -s "$d/refs/cfb.faa"  ] && label_refs "$d/refs/cfb.faa"  CFBCD  >> "$all"
        [ -s "$d/refs/bchl.faa" ] && label_refs "$d/refs/bchl.faa" BCHL   >> "$all"
        ;;
      mttB)
        echo "    fetching references ..."
        fetch "$d/refs/mttb.faa" 40 '"trimethylamine--corrinoid protein Co-methyltransferase"[Protein Name]'
        fetch "$d/refs/mtgb.faa" 40 '"glycine betaine--corrinoid protein Co-methyltransferase"[Protein Name] OR "MttB family protein"[Protein Name]'
        [ -s "$d/refs/mttb.faa" ] && label_refs "$d/refs/mttb.faa" TMA_MttB   >> "$all"
        [ -s "$d/refs/mtgb.faa" ] && label_refs "$d/refs/mtgb.faa" BETAINE_SF >> "$all"
        ;;
    esac
  fi

  local n; n=$(grep -c '^>' "$all")
  echo "    alignment input: $n sequences"
  if [ "$n" -lt 8 ]; then
    echo "    too few sequences for a meaningful tree — check NCBI connectivity"
    return 1
  fi

  # dereplicate identical headers, which break IQ-TREE
  run_tool seqkit rmdup -n "$all" -o "$d/${mk}_dedup.faa" 2>/dev/null || cp "$all" "$d/${mk}_dedup.faa"

  if [ ! -s "$d/${mk}_aln.faa" ]; then
    echo "    MAFFT ..."
    record_version "mafft" "${ENV_mafft:--}" mafft --version
    run_tool mafft --auto --thread "$THREADS" --quiet "$d/${mk}_dedup.faa" > "$d/${mk}_aln.faa" 2>/dev/null
  fi
  if [ ! -s "$d/${mk}_trim.faa" ]; then
    echo "    trimAl (gappyout) ..."
    run_tool trimal -in "$d/${mk}_aln.faa" -out "$d/${mk}_trim.faa" -gappyout 2>/dev/null \
      || cp "$d/${mk}_aln.faa" "$d/${mk}_trim.faa"
  fi

  if [ ! -s "$d/${mk}.treefile" ]; then
    echo "    IQ-TREE ..."
    run_tool iqtree3 -s "$d/${mk}_trim.faa" -m MFP -mset LG,WAG,JTT \
        -B 1000 -alrt 1000 -T AUTO --threads-max "$THREADS" \
        --prefix "$d/$mk" -redo > "$d/iqtree.log" 2>&1 \
      || echo "    !! IQ-TREE failed — see $d/iqtree.log"
  fi
  [ -s "$d/${mk}.treefile" ] && echo "    -> $d/${mk}.treefile"
}

build_marker mcrA "$V/mcrA_K00399_candidates.faa"
build_marker nifH "$V/nifH_K02588_candidates.faa"
build_marker mttB "$V/mttB_K14083_candidates.faa"

# =====================================================================
# VERDICTS — which reference group is each query sequence closest to?
# =====================================================================
python3 - "$M" "$WS2" <<'PY'
import sys, os, re, glob

M, WS2 = sys.argv[1], sys.argv[2]

def read_tree(p):
    return open(p).read().strip()

def neighbours(tree, leaf, window=400):
    """crude but transparent: the reference groups appearing nearest the query
    in the Newick string. A real assignment needs eyes on the tree — this is a
    triage aid, and it says so."""
    i = tree.find(leaf)
    if i < 0: return []
    seg = tree[max(0, i-window): i+window]
    return re.findall(r'(MCR|ACR|NIFH|CFBCD|BCHL|TMA_MttB|BETAINE_SF)\|', seg)

EXPECT = {
 'mcrA': ("MCR", "ACR",
          "MCR = true methanogen. ACR = alkane oxidiser; the claim would collapse."),
 'nifH': ("NIFH", "CFBCD",
          "NIFH = real nitrogenase. CFBCD = coenzyme F430 biosynthesis, NOT nitrogen fixation."),
 'mttB': ("TMA_MttB", "BETAINE_SF",
          "TMA_MttB = methylamine methanogenesis. BETAINE_SF = glycine betaine superfamily."),
}

for mk, (good, bad, note) in EXPECT.items():
    tf = os.path.join(M, mk, f'{mk}.treefile')
    print()
    print("=" * 88)
    print(f"  {mk} — TRIAGE")
    print("=" * 88)
    print(f"  {note}")
    if not os.path.exists(tf):
        print("  no tree built"); continue
    tree = read_tree(tf)
    queries = sorted(set(re.findall(r'QUERY\|([^,\):]+)', tree)))
    if not queries:
        print("  no query sequences found in the tree"); continue
    print()
    print(f"  {'query (MAG|contig|locus)':<52}{'nearest groups':<28}{'call':<10}")
    for q in queries:
        nb = neighbours(tree, 'QUERY|' + q)
        cnt = {g: nb.count(g) for g in set(nb)}
        top = sorted(cnt, key=lambda g: -cnt[g])
        call = 'CONFIRMED' if top and top[0] == good else ('REJECTED' if top and top[0] == bad else 'unclear')
        mag = q.split('|')[0]
        print(f"  {q[:52]:<52}{','.join(top[:3])[:28]:<28}{call:<10}")
    print()
    print("  This is Newick-adjacency triage, not a substitute for looking at the")
    print("  tree. Open the treefile in iTOL/FigTree and confirm the placement of")
    print("  every sequence before any of these calls goes into the thesis.")
PY

echo
echo "######## VALIDATION GATE ########"
for mk in mcrA nifH mttB; do
  check "$mk tree built" $([ -s "$M/$mk/${mk}.treefile" ] && echo 0 || echo 1)
  if [ -s "$M/$mk/${mk}_dedup.faa" ]; then
    nq=$(grep -c '^>QUERY' "$M/$mk/${mk}_dedup.faa" || true)
    nr=$(( $(grep -c '^>' "$M/$mk/${mk}_dedup.faa") - nq ))
    check "$mk has references (${nr}) alongside ${nq} queries" \
      $([ "$nr" -gt 5 ] && echo 0 || echo 1)
  fi
done

runlog "15_marker_trees complete"
echo
echo "######## 15 DONE ########"
echo "Open the trees to confirm —"
echo "the mcrA one decides whether 'methanogen' is a claim or a candidate."
