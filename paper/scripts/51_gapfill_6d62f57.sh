#!/usr/bin/env bash
# 51_gapfill.sh -- the second and final harvest pass.
#
#   bash scripts/51_gapfill.sh
#   -> workspace_2/REPORT_INPUTS_2.txt
#
# WHY A SECOND PASS. 50b did its job: it found every directory and showed me
# what is in them. But it protected itself two ways -- a ten-file dump cap per
# directory, and head+tail truncation on long tables -- and both fired exactly
# where the report needs whole tables:
#
#   * 43_smetana/run_43q/pairs and quads: 533 and 152 rows omitted. That is
#     THE community-modelling result of this chapter. Nothing can be said
#     about who feeds the methanogen without those rows in full.
#   * 29_transport, 23_extras, 34_unbinned, 35_kaiju, 36_gapclose,
#     41_env_checks, 13_profile SED: whole files never dumped at all because
#     the cap had already been reached in that directory.
#   * ~30 more tables truncated mid-way.
#
# So this script dumps a NAMED LIST of files IN FULL. No cap, no head+tail.
# Anything on the list that is not there prints [MISSING] -- which is a result,
# same as before. Nothing else is visited and nothing is written except the
# one output file.
#
# It is read-only apart from REPORT_INPUTS_2.txt, touches no environment, and
# runs no analysis.

set -uo pipefail

WS2="${WS2:-}"
if [ -z "$WS2" ]; then
  for c in ${WS2} \
           "$HOME/Documents/Metagenome_analysis_PhD/workspace_2"; do
    [ -d "$c" ] && { WS2="$c"; break; }
  done
fi
[ -n "$WS2" ] && [ -d "$WS2" ] || { echo "FATAL: set WS2=/path/to/workspace_2"; echo "STEP_STATUS: FAILED"; exit 1; }
OUT="$WS2/REPORT_INPUTS_2.txt"
HARDCAP=${HARDCAP:-600000}   # bytes; above this, head -400 only, and say so

exec > >(tee "$OUT") 2>&1

hr(){ printf '================================================================\n'; }
sec(){ printf '\n\n'; hr; printf '%s\n' "$*"; hr; }

full(){  # full RELATIVE_PATH
  local rel="$1" f="$WS2/$1" sz n
  if [ ! -e "$f" ]; then printf '\n>>> %s\n   [MISSING]\n' "$rel"; return; fi
  if [ ! -f "$f" ]; then printf '\n>>> %s\n   [NOT A FILE]\n' "$rel"; return; fi
  sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
  if [ "$sz" -eq 0 ]; then printf '\n>>> %s\n   [EMPTY]\n' "$rel"; return; fi
  if ! LC_ALL=C grep -Iq . "$f" 2>/dev/null; then
    printf '\n>>> %s  (%s bytes)\n   [BINARY] not dumped\n' "$rel" "$sz"; return; fi
  n=$(wc -l < "$f")
  if [ "$sz" -gt "$HARDCAP" ]; then
    printf '\n>>> %s   (%s bytes, %s lines) -- OVER HARD CAP, first 400 lines\n' "$rel" "$sz" "$n"
    head -n 400 "$f"
    printf '<<< TRUNCATED %s (over %s bytes)\n' "$rel" "$HARDCAP"
    return
  fi
  printf '\n>>> %s   (%s bytes, %s lines) -- FULL\n' "$rel" "$sz" "$n"
  cat "$f"
  printf '<<< end %s\n' "$rel"
}

hr
echo "ETOLIKO -- GAP-FILL HARVEST (51), the tables 50b truncated or capped"
echo "generated : $(date -Is)"
echo "ws2       : $WS2"
hr
echo
echo "Every file below is dumped IN FULL unless it is over $HARDCAP bytes."
echo "[MISSING] means the file is not there -- that is a result, not an error."

# ------------------------------------------------------------------ the result
sec "A. COMMUNITY MODELLING -- step 43q and its controls (THE headline)"
full 43_smetana/run_43q/pairs/pairs_detailed.tsv
full 43_smetana/run_43q/pairs/pairs_global.tsv
full 43_smetana/run_43q/pairs/communities.tsv
full 43_smetana/run_43q/quads/quads_detailed.tsv
full 43_smetana/run_43q/quads/quads_global.tsv
full 43_smetana/run_43q/quads/communities.tsv
full 43_smetana/run_43q/viable_pairs.json
full 43_smetana/run_43p/pairs/communities.tsv
full 43_smetana/run_43m/mcr_audit.json
full 43_smetana/run_43j/archaeon_fluxes.json
full 43_smetana/run_43i/currency_table.json
full 43_smetana/run_43i/currency_ids.json
full 43_smetana/run_43h/model_sets.json
full 43_smetana/run_43h/chosen.txt
full 43_smetana/run_43h/gapped/gap_global.tsv
full 43_smetana/run_43n/community/g_ETO_SED_global.tsv
full 43_smetana/run_43n/community/g_ETO_SED_FERM_global.tsv
full 43_smetana/run_43n/community/g_ETO_METH_global.tsv
full 43k_archaea/mtp4_gate.json

sec "B. THE MAG CATALOGUE AND ITS QUALITY, IN FULL"
full 12_results/MAG_CATALOGUE.md
full 12_results/MAG_PROFILES.md
full 06_binqc/MAG_QC.tsv
full 06_binqc/checkm2/quality_report.tsv
full 06_binqc/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv
full 06_binqc/rna_summary.tsv
full 06_binqc/union_membership.tsv
full 06_binqc/matched/ani.tsv
full 05_consensus/score_sweep.tsv
full 05_consensus/SED_DASTool_summary.tsv
full 23_extras/genome_stats.tsv
full 23_extras/rrna_copy_number.tsv

sec "C. TAXONOMY AND PHYLOGENY"
full 08_taxonomy/gtdbtk/gtdbtk.json
full 12_phylogeny/ani/ani_all.tsv
full 41_env_checks/mag_x_mag_ani.tsv
full 41_env_checks/in_situ_10.txt
full 40_integration/skani_raw.tsv
full 40_integration/INTEGRATION.txt

sec "D. FUNCTION -- the full per-genome module tables"
full 10_function/pathway_summary.tsv
full 10_function/syntrophy_handoff.tsv
full 11_kofam/ko_modules.tsv
full 11_kofam/ko_targets.tsv
full 45_kofam_full/b12_genomewide.json
full 44_loose_ends/b12_kofam.json
full 44_loose_ends/kofam_inventory.json
full 48_sulfur/sulfur_kos.json

sec "E. AUXOTROPHY AND THE INTERACTION NETWORK, IN FULL"
full 19_auxotrophy/auxotrophies.tsv
full 19b_auxotrophy/auxotrophies.tsv
full 19b_auxotrophy/routes.tsv
full 19b_auxotrophy/marker_audit.tsv
full 21_network/crossfeeding.tsv
full 21_network/competition_matrix.tsv
full 29_transport/edge_transport.tsv
full 29_transport/transport_summary.tsv
full 36_gapclose/transporter_matrix.tsv
full 36_gapclose/transporter_kos.tsv
full 36_gapclose/kos_all.tsv
full 36_gapclose/VERDICT.txt

sec "F. THE VERDICT TEXTS -- each step's own conclusion, in full"
full 14_comparative/BIOGEOGRAPHY.txt
full 14_comparative/summary.tsv
full 14_comparative/recruitment_matrix.tsv
full 16_environment/ENVIRONMENT.txt
full 34_unbinned/VERDICT.txt
full 34_unbinned/LENGTH_BIAS.txt
full 34b_methylamine/VERDICT.txt
full 35_kaiju/VERDICT.txt
full 35_kaiju/ASSIGNMENT_RANK.txt
full 36_gapclose/CARGO.txt
full 36_gapclose/CAZY_CROSS.txt
full 36_gapclose/ecf/ECF_VERDICT.txt
full 39_loose_ends/LOOSE_ENDS.txt
full 39_loose_ends/SECRETOME_BOUND.txt
full 39_loose_ends/provirus/REASSIGN.txt
full 32_strain/ANATOMY.txt
full 32_strain/REPEAT_CONTROL.txt
full 31_readlevel/VERDICT.txt
full 13_profile/silva/verify/VERDICT.txt
full 13_profile/silva/verify/RESOLUTION.txt
full 13_profile/silva/verify/DIAGNOSIS.txt

sec "G. COMMUNITY COMPOSITION -- the sediment side, which 50b's cap hid"
full 13_profile/kraken/SED_bracken_P.tsv
full 13_profile/kraken/SED_bracken_C.tsv
full 13_profile/kraken/SED_bracken_O.tsv
full 13_profile/kraken/SED_bracken_F.tsv
full 13_profile/kraken/SED_bracken_G.tsv
full 13_profile/kraken/SED_bracken_S.tsv
full 13_profile/SED_assembly_composition.tsv
full 13_profile/LIQ_assembly_composition.tsv
full 13_profile/silva/lineages_not_in_MAGs.tsv
full 22_singlem/SED_profile.tsv
full 22_singlem/LIQ_profile.tsv
full 24_deep_taxonomy/singlem_ranks.tsv
full 24_deep_taxonomy/silva_len400/lineages_not_in_MAGs.tsv

sec "H. SUBSTRATES, CAZYmes, SECRETION, VIRUSES"
full 23_extras/pul/pul_summary.tsv
full 23_extras/pul/secreted_vs_pul.tsv
full 23_extras/secretome.tsv
full 23_extras/secretion/secretion_systems.tsv
full 23_extras/secretion/motif_substrates.tsv
full 23_extras/t9ss/t9ss_machinery.tsv
full 23_extras/t9ss/cazyme_endo_exo.tsv
full 20_substrates/substrate_matrix.tsv
full 20_substrates/SUBSTRATES.txt
full 20_substrates/ko_targets.tsv
full 18_viruses/mags_genomad/all_mags_contigs_summary/all_mags_contigs_virus_summary.tsv
full 18_viruses/provirus_host_links.tsv
full 39_loose_ends/provirus/mag_proviruses.tsv

sec "I. ANYTHING NAMED ABOVE THAT WAS NOT FOUND -- and a last sweep"
echo "Files whose names I guessed and that do not exist show [MISSING] above."
echo
echo "Sweep: every .txt/.tsv in 20_substrates and 29_transport (both were capped):"
for d in 20_substrates 29_transport 23_extras/secretion 36_gapclose; do
  echo "  -- $d --"
  find "$WS2/$d" -maxdepth 2 -type f \( -name '*.tsv' -o -name '*.txt' \) 2>/dev/null |
    sort | while read -r f; do printf '   %10s  %s\n' "$(stat -c %s "$f")" "${f#$WS2/}"; done
done

hr
echo "output : $OUT"
echo "size   : $(du -h "$OUT" 2>/dev/null | cut -f1)"
echo
echo "If REPORT_INPUTS_2.txt is too large to handle in one piece:"
echo "   split -b 4m $OUT ${OUT}.part_"
hr
echo "STEP_STATUS: COMPLETE  $(date -Is)"
