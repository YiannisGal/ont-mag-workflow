#!/usr/bin/env bash
# 50b_harvest_complete.sh -- dump the primary output of EVERY directory in
# workspace_2 into one file, so the findings report is written from
# the data.
#
#   bash scripts/50b_harvest_complete.sh
#   -> workspace_2/REPORT_INPUTS.txt
#
# WHY THIS EXISTS. Version 1 of the harvest covered about twenty directories.
# The real listing has sixty. Anything not harvested would have had to be
# written from memory -- and eight separate places have already been caught
# where a remembered number was wrong (Methanosarcina 3 vs 2, "0 B12
# consumers", "no mcr in vamb_017", "archaeal genomes: 0", the E-value read as a
# score, and so on). So this walks the directory list literally, in the order it
# was given, and prints MISSING where a directory or a file is not there. A
# MISSING line is a real result: it means that analysis was never run, or was
# run somewhere else, and the report must not claim it.
#
# WHAT IT DOES NOT DO. It never writes into workspace_1, never touches a conda
# env, never deletes anything, and never runs an analysis. It is read-only
# except for the one output file. Anomalous directories (the garbage-named ones)
# are INSPECTED and REPORTED, not removed -- you decide.
#
# SIZE. Big tables are truncated head+tail with the omitted-line count shown, so
# the file stays manageable. Anything truncated is named, so the
# full version of any single table can be pulled afterwards.

set -uo pipefail

# ---------------------------------------------------------------- locate ws2
WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
if [ -z "$WS2" ] || [ ! -d "$WS2" ]; then
  echo "FATAL: cannot find workspace_2. Run with:  WS2=/path/to/workspace_2 bash $0"
  echo "STEP_STATUS: FAILED"; exit 1
fi
OUT="$WS2/REPORT_INPUTS.txt"

MAXFULL=${MAXFULL:-70}     # tables at or under this many lines are dumped whole
MAXBYTES=${MAXBYTES:-600000}   # never cat a text file larger than this
DUMPCAP=${DUMPCAP:-10}     # generic text files dumped per directory, after key files

# ---------------------------------------------------------------- helpers
hr(){ printf '================================================================\n'; }
sec(){ printf '\n\n'; hr; printf '%s\n' "$*"; hr; }
sub(){ printf '\n---------------- %s\n' "$*"; }

is_text(){ [ -f "$1" ] && LC_ALL=C grep -Iq . "$1" 2>/dev/null; }

fsize(){ stat -c %s "$1" 2>/dev/null || echo 0; }

dump(){   # dump FILE [maxlines]
  local f="$1" mx="${2:-$MAXFULL}" n sz
  if [ ! -e "$f" ]; then printf '   [MISSING]  %s\n' "$f"; return 1; fi
  if [ ! -f "$f" ]; then printf '   [NOT A FILE] %s\n' "$f"; return 1; fi
  sz=$(fsize "$f")
  if [ "$sz" -eq 0 ]; then printf '   [EMPTY]    %s\n' "$f"; return 0; fi
  if ! is_text "$f"; then printf '   [BINARY]   %s  (%s bytes) -- not dumped\n' "$f" "$sz"; return 0; fi
  if [ "$sz" -gt "$MAXBYTES" ]; then
    n=$(wc -l < "$f")
    printf '\n>>> %s   (%s bytes, %s lines) -- TOO LARGE, head only\n' "$f" "$sz" "$n"
    head -n 15 "$f"
    printf '<<< truncated (%s bytes). Ask for this one separately if needed.\n' "$sz"
    return 0
  fi
  n=$(wc -l < "$f")
  printf '\n>>> %s   (%s bytes, %s lines)\n' "$f" "$sz" "$n"
  if [ "$n" -le "$mx" ]; then
    cat "$f"
  else
    head -n $((mx - 8)) "$f"
    printf '   ...... [%s lines omitted] ......\n' "$(( n - (mx - 8) - 6 ))"
    tail -n 6 "$f"
  fi
  printf '<<< end %s\n' "$(basename "$f")"
}

# key files per directory: always dumped, found anywhere beneath the directory,
# and reported MISSING by name if absent. Lesson from step 44: the file you want
# is often one level down (11_kofam/hits/), so never glob only the top level.
keyfiles(){
  case "$1" in
    00_input)         echo "manifest.tsv samples.tsv README.txt" ;;
    01_qc)            echo "nanoplot_summary.tsv qc_summary.tsv NanoStats.txt summary.tsv" ;;
    02_mapping)       echo "mapping_stats.tsv flagstat.txt coverage.tsv depth_summary.tsv" ;;
    03_assembly)      echo "assembly_stats.tsv assembly_summary.tsv quast_report.tsv report.tsv contig_stats.tsv" ;;
    03_binning|03_bins) echo "binning_summary.tsv bins_summary.tsv" ;;
    04_polish)        echo "polish_summary.tsv medaka_summary.tsv" ;;
    05_binqc|06_binqc) echo "checkm2_quality_report.tsv quality_report.tsv checkm_summary.tsv" ;;
    05_consensus)     echo "consensus_summary.tsv dastool_summary.tsv DASTool_summary.tsv" ;;
    06_taxonomy|08_taxonomy|24_deep_taxonomy) echo "gtdbtk.bac120.summary.tsv gtdbtk.ar53.summary.tsv taxonomy_summary.tsv" ;;
    07_mags)          echo "mag_table.tsv mags.tsv naming_map.tsv mag_inventory.tsv" ;;
    08_annotation|09_annotation) echo "annotation_summary.tsv prokka_summary.tsv bakta_summary.tsv" ;;
    09_coverage)      echo "coverage_summary.tsv breadth.tsv depth.tsv" ;;
    10_function)      echo "function_summary.tsv kegg_summary.tsv pathways.tsv" ;;
    10_phylogenomics|12_phylogeny|15_marker_trees) echo "tree_summary.tsv concat.treefile iqtree.log" ;;
    11_kofam)         echo "kofam_summary.tsv ko_table.tsv" ;;
    11_markers)       echo "marker_summary.tsv" ;;
    12_results)       echo "results_summary.tsv summary.tsv" ;;
    13_profile)       echo "profile_summary.tsv abundance.tsv" ;;
    14_comparative|17_compgen) echo "comparative_summary.tsv ani_summary.tsv fastani.tsv" ;;
    16_environment)   echo "environment_summary.tsv medium.tsv" ;;
    18_viruses)       echo "virus_summary.tsv checkv_quality_summary.tsv" ;;
    19_auxotrophy|19b_auxotrophy) echo "auxotrophy_summary.tsv auxotrophies.tsv" ;;
    20_substrates)    echo "substrate_summary.tsv substrates.tsv" ;;
    21_network)       echo "network_summary.tsv interactions.tsv" ;;
    22_singlem)       echo "singlem_summary.tsv otu_table.tsv profile.tsv" ;;
    26_crossmap)      echo "cross_abundance.tsv crossmap_summary.tsv" ;;
    27_rrna_corrected) echo "rrna_summary.tsv barrnap_summary.tsv" ;;
    28_crispr)        echo "crispr_summary.tsv" ;;
    29_transport)     echo "transport_summary.tsv" ;;
    30_dbcan)         echo "dbcan_summary.tsv cazy_summary.tsv" ;;
    31_readlevel)     echo "readlevel_summary.tsv" ;;
    32_strain)        echo "strain_summary.tsv" ;;
    33_antismash)     echo "antismash_summary.tsv bgc_summary.tsv" ;;
    34_unbinned)      echo "unbinned_summary.tsv" ;;
    34b_methylamine)  echo "methylamine_summary.tsv" ;;
    35_kaiju)         echo "kaiju_summary.tsv" ;;
    36_gapclose)      echo "gapclose_summary.tsv" ;;
    38_gems|42_gems)  echo "gems_summary.tsv models_summary.tsv" ;;
    38_preflight|41_env_checks) echo "preflight.tsv env_checks.tsv" ;;
    39_loose_ends)    echo "loose_ends.tsv" ;;
    40_integration)   echo "integration_summary.tsv" ;;
    43_smetana)       echo "smetana_summary.tsv pc_detailed.tsv" ;;
    43k_archaea)      echo "mtp4_gate.txt" ;;
    44_loose_ends)    echo "in_situ_ten.tsv ani_pairs.tsv taxa_counts.tsv" ;;
    45_kofam_full)    echo "ko_per_genome.tsv kofam_hits_above_threshold.tsv b12_summary.tsv" ;;
    47_master)        echo "MASTER.tsv" ;;
    48_sulfur)        echo "sulfur_summary.tsv dsr_summary.tsv" ;;
    49_dsr)           echo "dsr_contigs.tsv dsr_contig_proteins.faa depth_control.tsv" ;;
    50_figures)       echo "" ;;
    *)                echo "" ;;
  esac
}

# ---------------------------------------------------------------- start
exec > >(tee "$OUT") 2>&1

hr
echo "ETOLIKO LAGOON -- workspace_2 COMPLETE HARVEST (50b)"
echo "generated : $(date -Is)"
echo "host      : $(hostname)"
echo "user      : $(whoami)"
echo "ws2       : $WS2"
echo "real path : $(readlink -f "$WS2")"
hr
echo
echo "HOW TO READ THIS FILE"
echo "  [MISSING]  the file named is not there -> that analysis did not produce it"
echo "  [EMPTY]    the file exists with zero bytes -> the step ran and produced nothing"
echo "  [BINARY]   not dumped, size reported"
echo "  >>> ... <<<  a dumped table; if truncated the omitted line count is shown"
echo
echo "Every directory from the ls you sent is visited below, in that order."
echo "Nothing here is interpreted. Interpretation happens in the report."

# ---------------------------------------------------------------- part 0
sec "PART 0 -- TOP LEVEL"

sub "0.1  ls -la"
ls -la "$WS2"

sub "0.2  disk"
df -h "$WS2" 2>/dev/null
echo
echo "du -sh of workspace_2 (may take a moment):"
du -sh "$WS2" 2>/dev/null

sub "0.3  every directory with its size and file count"
printf '%-34s %8s %8s\n' DIRECTORY SIZE FILES
for d in "$WS2"/*/; do
  [ -d "$d" ] || continue
  b=$(basename "$d")
  s=$(du -sh "$d" 2>/dev/null | cut -f1)
  n=$(find "$d" -type f 2>/dev/null | wc -l)
  printf '%-34s %8s %8s\n' "$b" "$s" "$n"
done

sub "0.4  config"
for f in 00_config.sh 00_config_patch.sh env_resolved.sh versions.tsv; do
  dump "$WS2/$f" 120
done

sub "0.5  run logs"
for f in runlog.md RUNLOG.md; do
  if [ -f "$WS2/$f" ]; then
    echo "(both runlog.md and RUNLOG.md exist -- same file on a case-insensitive mount, or two files?)"
    dump "$WS2/$f" 200
  else
    printf '   [MISSING]  %s\n' "$WS2/$f"
  fi
done
echo
echo "inode check (identical inode = one file seen twice):"
stat -c '%i  %n' "$WS2/runlog.md" "$WS2/RUNLOG.md" 2>&1

sub "0.6  labnotes/"
if [ -d "$WS2/labnotes" ]; then
  ls -la "$WS2/labnotes"
else
  echo "   [MISSING]  $WS2/labnotes"
fi

sub "0.7  scripts/ (names and md5, so we can tell which version ran)"
if [ -d "$WS2/scripts" ]; then
  find "$WS2/scripts" -maxdepth 1 -type f | sort | while read -r f; do
    printf '%s  %10s  %s\n' "$(md5sum "$f" | cut -c1-32)" "$(fsize "$f")" "$(basename "$f")"
  done
else
  echo "   [MISSING]  $WS2/scripts"
fi

sub "0.8  logs/ (names, sizes, last line of each)"
if [ -d "$WS2/logs" ]; then
  find "$WS2/logs" -maxdepth 2 -type f | sort | while read -r f; do
    printf '\n%10s  %s\n' "$(fsize "$f")" "$f"
    if is_text "$f"; then
      echo "        last: $(tail -n 1 "$f" 2>/dev/null | cut -c1-160)"
    fi
  done
else
  echo "   [MISSING]  $WS2/logs"
fi

# ---------------------------------------------------------------- part 1
sec "PART 1 -- EVERY DIRECTORY, IN ORDER"

DIRS="
00_input
00_smoketest
01_qc
02_mapping
03_assembly
03_binning
04_polish
05_binqc
05_consensus
06_binqc
06_taxonomy
07_mags
08_annotation
08_taxonomy
09_annotation
09_coverage
10_function
10_phylogenomics
11_kofam
11_markers
12_phylogeny
12_results
13_profile
14_comparative
15_marker_trees
16_environment
17_compgen
18_viruses
19_auxotrophy
19b_auxotrophy
20_substrates
21_network
22_singlem
23_extras
24_deep_taxonomy
26_crossmap
27_rrna_corrected
28_crispr
29_transport
30_dbcan
31_readlevel
32_strain
33_antismash
34_unbinned
34b_methylamine
35_kaiju
36_gapclose
38_gems
38_preflight
39_loose_ends
40_integration
41_env_checks
42_gems
43k_archaea
43_smetana
44_loose_ends
45_kofam_full
47_master
48_sulfur
49_dsr
50_figures
bin
envs
tmp
"

MISSING_DIRS=""
for name in $DIRS; do
  D="$WS2/$name"
  sec "1.$name"
  if [ ! -d "$D" ]; then
    echo "   [MISSING DIRECTORY]  $D"
    echo "   -> nothing from this step exists; the report will say so."
    MISSING_DIRS="$MISSING_DIRS $name"
    continue
  fi

  echo "path : $D"
  echo "size : $(du -sh "$D" 2>/dev/null | cut -f1)"
  echo "files: $(find "$D" -type f 2>/dev/null | wc -l)   dirs: $(find "$D" -mindepth 1 -type d 2>/dev/null | wc -l)"
  echo "mtime: $(stat -c %y "$D" 2>/dev/null)"

  # envs/ and bin/ and tmp/: inventory only, never dump
  case "$name" in
    envs|bin|tmp)
      sub "inventory only (no dumps)"
      find "$D" -maxdepth 2 2>/dev/null | sort | head -60
      t=$(find "$D" 2>/dev/null | wc -l)
      [ "$t" -gt 60 ] && echo "   ... $((t-60)) more entries"
      continue ;;
  esac

  sub "contents (depth 2)"
  find "$D" -maxdepth 2 -mindepth 1 -printf '%y %10s  %P\n' 2>/dev/null | sort -k3 | head -100
  t=$(find "$D" -maxdepth 2 -mindepth 1 2>/dev/null | wc -l)
  [ "$t" -gt 100 ] && echo "   ... $((t-100)) more entries at depth<=2"

  # census by extension, at ANY depth, so a file type I did not think to look for
  # (.tbl, .gff, .treefile, .gbk) is still visible in the count even when nothing
  # matching it gets dumped. Step 44 missed 200 .tbl files this way.
  sub "file types at any depth"
  find "$D" -type f -printf '%f\n' 2>/dev/null |
    awk '{n=split($0,a,"."); e=(n>1? a[n] : "(no extension)"); c[e]++}
         END{for(k in c) printf "   %8d  .%s\n", c[k], k}' | sort -k1 -nr | head -25

  # ---- key files, searched at ANY depth
  KF="$(keyfiles "$name")"
  if [ -n "$KF" ]; then
    sub "expected primary outputs"
    for k in $KF; do
      hits=$(find "$D" -type f -name "$k" 2>/dev/null | sort)
      if [ -z "$hits" ]; then
        printf '   [MISSING]  %s  (searched all of %s)\n' "$k" "$name"
      else
        echo "$hits" | while read -r h; do dump "$h"; done
      fi
    done
  fi

  # ---- generic sweep: other small text tables not already dumped
  sub "other text output in this directory"
  shown=0
  find "$D" -type f \
       \( -name '*.tsv' -o -name '*.csv' -o -name '*.txt' -o -name '*.tab' \
          -o -name '*.json' -o -name '*.md' -o -name '*summary*' -o -name '*report*' \) \
       -size -400k 2>/dev/null | sort | while read -r f; do
    base=$(basename "$f")
    skip=0
    for k in $KF; do [ "$base" = "$k" ] && skip=1; done
    [ "$skip" -eq 1 ] && continue
    shown=$((shown+1))
    if [ "$shown" -gt "$DUMPCAP" ]; then
      printf '   (dump cap %s reached in %s -- remaining text files listed only)\n' "$DUMPCAP" "$name"
      break
    fi
    dump "$f" 45
  done
  # list everything text-ish so nothing is invisible even when capped
  echo
  echo "all text-ish files in $name (name, size, lines):"
  find "$D" -type f \( -name '*.tsv' -o -name '*.csv' -o -name '*.txt' -o -name '*.tab' \
       -o -name '*.json' -o -name '*.md' -o -name '*.log' \) 2>/dev/null | sort | while read -r f; do
    if is_text "$f"; then
      printf '   %10s  %8s lines  %s\n' "$(fsize "$f")" "$(wc -l < "$f")" "${f#$WS2/}"
    else
      printf '   %10s  %8s        %s\n' "$(fsize "$f")" "-" "${f#$WS2/}"
    fi
  done | head -120

  # ---- sequence inventory, counted not dumped
  nfa=$(find "$D" -type f \( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' -o -name '*.faa' \) 2>/dev/null | wc -l)
  if [ "$nfa" -gt 0 ]; then
    sub "sequence files ($nfa) -- counted, not dumped"
    find "$D" -type f \( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' -o -name '*.faa' \) 2>/dev/null |
      sort | head -40 | while read -r f; do
        printf '   %10s  %7s seqs  %s\n' "$(fsize "$f")" "$(grep -c '^>' "$f" 2>/dev/null || echo '?')" "${f#$WS2/}"
      done
    [ "$nfa" -gt 40 ] && echo "   ... $((nfa-40)) more"
  fi
done

# ---------------------------------------------------------------- part 2
sec "PART 2 -- ANOMALIES IN THE TOP LEVEL"
echo "These are not analysis directories. Most look like a shell variable that"
echo "expanded to nothing and a mkdir that then made a directory named after"
echo "whatever junk followed. They are INSPECTED here, not removed."
echo

ANOM="78CC AAIEEEDHMLHRRQ AB BB BIHCJIKKLKLDDDD CECCBCHHHHFGIJLMHE DEGHJJLJJNNEEFHNAIIJMF DJNLGEEEHBEA GFGGEECDDEHHJHJJML GIIGFE IJJIJDBBABBCBBDAACBC 5 91128.pts-2.MGeMServer _quarantine_20260825"
for a in $ANOM; do
  P="$WS2/$a"
  sub "$a"
  if [ ! -e "$P" ]; then echo "   [ABSENT NOW]"; continue; fi
  echo "   type : $(stat -c %F "$P" 2>/dev/null)"
  echo "   size : $(du -sh "$P" 2>/dev/null | cut -f1)"
  echo "   mtime: $(stat -c %y "$P" 2>/dev/null)"
  if [ -d "$P" ]; then
    n=$(find "$P" -type f 2>/dev/null | wc -l)
    echo "   files: $n"
    find "$P" -maxdepth 2 2>/dev/null | sort | head -25
  elif [ -f "$P" ]; then
    if is_text "$P"; then head -n 10 "$P"; else echo "   (binary)"; fi
  fi
done

sub "stray files in the top level"
for f in scripts43p_pairwise_78fe4f9.sh 43p_pairwise.sh ctd.hmm hmm_PGAP.HMM.tgz 00_inventory.sh 00_config.sh.bak_20260823_1257 00_config.sh.bak_20260823_1745; do
  P="$WS2/$f"
  if [ -e "$P" ]; then
    printf '   %10s  %s  %s\n' "$(fsize "$P")" "$(stat -c %y "$P" | cut -c1-16)" "$f"
  else
    printf '   [ABSENT]   %s\n' "$f"
  fi
done
echo
echo "scripts43p_pairwise_78fe4f9.sh is an scp with a missing slash"
echo "(scp ... server:${WS2}/scripts43p... instead of .../scripts/43p...)."
echo "If it is byte-identical to scripts/43p_pairwise_78fe4f9.sh it is safe to remove:"
if [ -f "$WS2/scripts43p_pairwise_78fe4f9.sh" ]; then
  md5sum "$WS2/scripts43p_pairwise_78fe4f9.sh" 2>/dev/null
  find "$WS2/scripts" -maxdepth 1 -name '43p*' -exec md5sum {} \; 2>/dev/null
fi

sub "hmm_PGAP / ctd.hmm -- where did these come from"
for p in "$WS2/hmm_PGAP" "$WS2/ctd.hmm" "$WS2/hmm_PGAP.HMM.tgz"; do
  [ -e "$p" ] || { echo "   [ABSENT] $p"; continue; }
  echo "   $p  $(du -sh "$p" 2>/dev/null | cut -f1)  $(stat -c %y "$p" | cut -c1-16)"
done

# ---------------------------------------------------------------- part 3
sec "PART 3 -- NUMBERS I AM CARRYING IN MEMORY, RE-DERIVED FROM THE FILES"
echo "Each block below re-computes a number the report would otherwise state from"
echo "memory. If a block prints nothing or MISSING, the report will not make that"
echo "claim. Do not fix anything here -- the point is to see what is actually true."

M="$(find "$WS2/47_master" -type f -name MASTER.tsv 2>/dev/null | head -1)"
sub "3.1  MASTER.tsv columns"
if [ -n "$M" ]; then
  head -1 "$M" | tr '\t' '\n' | nl
else
  echo "   [MISSING] MASTER.tsv"
fi

sub "3.2  genome counts by domain"
if [ -n "$M" ]; then
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;
                    if(!("domain" in h)){print "   [NO domain COLUMN]"; bad=1} next}
              bad{next}
              {c[$h["domain"]]++} END{for(k in c) printf "   %-12s %s\n",k,c[k]}' "$M"
  echo "   total rows: $(( $(wc -l < "$M") - 1 ))"
else echo "   [MISSING]"; fi

sub "3.3  every genome by sed_breadth, with its in_situ flag"
if [ -n "$M" ]; then
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
    function g(k){ return (k in h)? $h[k] : "NA" }
    {printf "%s\t%s\t%s\t%s\t%s\n",g("mag_id"),g("sed_breadth"),g("in_situ_ten"),g("domain"),g("genus")}
    ' "$M" | sort -t$'\t' -k2,2gr | nl
  echo "   (columns: mag_id, sed_breadth, in_situ_ten, domain, genus; NA = column absent)"
else echo "   [MISSING]"; fi

sub "3.4  genus counts (exact field match, not substring)"
if [ -n "$M" ]; then
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;
                    if(!("genus" in h)){print "   [NO genus COLUMN]"; bad=1} next}
              bad{next}
              {g=$h["genus"]; if(g!=""&&g!="-"&&g!="NA")c[g]++}
              END{for(k in c) printf "   %-30s %s\n",k,c[k]}' "$M" | sort -k2 -nr
else echo "   [MISSING]"; fi

sub "3.5  SMETANA result table (43q / 43_smetana)"
for f in $(find "$WS2/43_smetana" "$WS2/43p_pairwise" -type f \( -name '*detailed*.tsv' -o -name '*smetana*.tsv' -o -name '*pairwise*.tsv' \) 2>/dev/null | sort); do
  dump "$f" 80
done

sub "3.6  KOfam totals"
for f in $(find "$WS2/45_kofam_full" -maxdepth 2 -type f -name '*.tsv' -size -400k 2>/dev/null | sort); do
  dump "$f" 45
done

sub "3.7  dsrAB evidence (48_sulfur, 49_dsr)"
for f in $(find "$WS2/48_sulfur" "$WS2/49_dsr" -type f \( -name '*.tsv' -o -name '*.txt' \) -size -400k 2>/dev/null | sort); do
  dump "$f" 60
done

sub "3.8  cross-mapping abundance (source of sed_breadth)"
for f in $(find "$WS2/26_crossmap" -type f -name '*.tsv' -size -400k 2>/dev/null | sort); do
  dump "$f" 70
done

sub "3.9  GTDB-Tk summaries, both domains"
for f in $(find "$WS2" -type f \( -name 'gtdbtk.bac120.summary.tsv' -o -name 'gtdbtk.ar53.summary.tsv' \) 2>/dev/null | sort); do
  dump "$f" 45
done

sub "3.10  CheckM2 quality, wherever it lives"
for f in $(find "$WS2" -type f \( -name 'quality_report.tsv' -o -name 'checkm2*.tsv' \) 2>/dev/null | sort); do
  dump "$f" 45
done

# ---------------------------------------------------------------- part 4
sec "PART 4 -- MANIFEST"
echo "Directories in the list that do not exist:"
if [ -z "$MISSING_DIRS" ]; then echo "   (none -- every directory in the ls is present)"
else for d in $MISSING_DIRS; do echo "   $d"; done; fi
echo
echo "Directories present in workspace_2 but NOT in my list (so not harvested above):"
# $DIRS is newline-separated; unquoted echo collapses it to single spaces so the
# case pattern can match. Quoting it here would leave the newlines in and make
# every harvested directory look unharvested -- which is exactly what v1 of this
# block did.
DIRS_FLAT="$(echo $DIRS)"
KNOWN=" $DIRS_FLAT $ANOM labnotes logs scripts "
unlisted=0
for d in "$WS2"/*/; do
  b=$(basename "$d")
  case "$KNOWN" in *" $b "*) ;; *) echo "   $b"; unlisted=$((unlisted+1)) ;; esac
done
[ "$unlisted" -eq 0 ] && echo "   (none -- every directory present was harvested or inspected)"
echo
echo "output file : $OUT"
echo "output size : $(du -h "$OUT" 2>/dev/null | cut -f1)"
echo
echo "If REPORT_INPUTS.txt is too large to handle in one piece, run:"
echo "   split -b 4m $OUT ${OUT}.part_"
hr
echo "STEP_STATUS: COMPLETE  $(date -Is)"
