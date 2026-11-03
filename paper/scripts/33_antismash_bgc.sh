#!/usr/bin/env bash
# =============================================================================
# 33_antismash_bgc.sh  --  biosynthetic gene clusters: what do these genomes
#                          make that is not primary metabolism?
#
# NUMBERING
#   The roadmap called this step 31. Step 31 was spent on the read-level
#   methylotrophy test, which arose from step 30 and was never on the roadmap.
#   That drift is recorded in remaining_analyses_roadmap.md under Rule 16.
#   antiSMASH is 33. Nothing that already ran was renamed.
#
# THE QUESTION
#   Everything in this chapter so far is metabolism the organism does to feed
#   itself: CAZymes, transporters, methyltransferases, sulfate reduction.
#   Secondary metabolism is the other layer — what an organism makes to act on
#   the organisms around it. In a community held at 29 genomes with a
#   cross-feeding network (step 21) and a virus/defence landscape (step 28),
#   bacteriocins and RiPPs are the competition arm of that same story.
#
#   The honest expectation is a MODEST yield. Anaerobic sediment bacteria are
#   not Streptomyces. But a genuine near-zero in a community this novel is
#   itself a reportable result — provided it is a TESTED zero and not a
#   not-looked-for zero, which is the distinction step 20b and step 31 were
#   built around.
#
# TWO PATHS, AND THE SECOND IS NOT THE FIRST
#   PATH A  antiSMASH, if it is installed AND its databases are present.
#           This is the real analysis: cluster boundaries, rule-based product
#           classes, MIBiG comparison.
#   PATH B  a Pfam core-domain survey using the HMM library Bakta already
#           ships. This is a FLOOR, not a substitute. It finds the biosynthetic
#           machinery; it cannot draw a cluster boundary, cannot name a
#           product, and will miss every RiPP class whose only signature is a
#           short unannotated precursor peptide. Labelled as such everywhere it
#           is reported.
#
# WHY THE FALLBACK NEEDS A DISCRIMINATOR, NOT A DOMAIN COUNT
#   Counting core domains naively gives a wrong answer with high confidence:
#     - AMP-binding is a huge superfamily. Every fatty-acid-CoA ligase has one.
#     - PP-binding (ACP) is in fatty acid synthesis in every bacterium alive.
#     - Ketoacyl-synthase is FabB and FabF — again, primary metabolism.
#   The features that actually separate secondary from primary metabolism are:
#     1. the CONDENSATION domain, which is essentially NRPS-specific;
#     2. ARCHITECTURE — bacterial FAS-II is dissociated (one domain per protein,
#        ~400 aa), while type-I PKS and NRPS are megasynthases carrying several
#        core domains in ONE ORF of 1000-5000 aa.
#   So PKS is called on multidomain-in-one-ORF plus length, never on a KS hit.
#   The length distribution of KS-carrying proteins is printed so the reader can
#   see the discriminator separating two populations rather than trust it.
#
# MODEL SELECTION IS DERIVED, NOT ASSERTED  (rule 1)
#   No Pfam accession is written into this script from memory. Models are
#   selected by matching their own DESC line against keyword patterns, and every
#   model that matched is printed so the selection can be audited. If the
#   library is pressed-only, DESC is recovered from hmmscan's target-description
#   column instead — the same text, from the other direction.
#
# THE CONTROLS
#   POSITIVE  PP-binding and AMP-binding must be found in nearly every BACTERIAL
#             genome, because fatty acid synthesis is not optional. If they are
#             not, the search is broken and NO VERDICT is printed. This control
#             comes from the data, like step 31's archaeal ribosomal control.
#   ARCHAEA   are scored SEPARATELY and are NOT part of that control. Archaea
#             build isoprenoid ether lipids and do not run bacterial FAS-II, so
#             an ACP/KS absence in the eight archaeal genomes is EXPECTED and is
#             not evidence of a failed search. Applying a bacterial control set
#             to archaea is the exact error of steps 19b and 23e.
#   NEGATIVE  the methanogens should carry few or no NRPS condensation domains.
#             A methanogen full of NRPS is a sign the filter is too loose.
#
# IN SITU WEIGHTING
#   Step 26 showed six culture genomes are real sediment residents and the rest
#   are enrichment-specific. A BGC in a genome at 0.02 sediment breadth says
#   nothing about the lagoon. Counts are therefore split by sediment breadth and
#   the headline number is the in-situ one.
#
# INPUT   09_annotation/<genome>/<genome>.{gbff,faa,tsv}   07_mags/MAG_TAXONOMY.tsv
#         26_crossmap/cross_abundance.tsv                  (optional, for weighting)
# OUTPUT  33_antismash/bgc_regions.tsv | core_domains.tsv | BGC_SUMMARY.tsv
#         33_antismash/VERDICT.txt
# RUNTIME antiSMASH 3-10 min/genome. Pfam fallback 20-90 min, cached.
# USAGE   ws2 && bash scripts/33_antismash_bgc.sh 2>&1 | tee logs/33_antismash.log
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/../00_config.sh"
export STEP="33_antismash"
export TMPDIR="${TMPDIR:-${TMPDIR}}"; mkdir -p "$TMPDIR"

A="$WS2/09_annotation"; F="$WS2/07_mags"; X="$WS2/26_crossmap"
B="$WS2/33_antismash"; mkdir -p "$B" "$B/as" "$B/hmm"

# ---------------------------------------------------------------- PREFLIGHT
# The first two attempts at this step left output behind: a BGC_SUMMARY.tsv and
# a VERDICT.txt from the Pfam run whose control failed, and 29 antiSMASH output
# directories holding nothing but the log file that caused the abort. None of it
# is a measurement. Derived files are therefore removed at the START of every
# run, not merely overwritten at the end — if this run dies halfway, a stale
# summary from a previous attempt must not be sitting there looking current.
rm -f "$B/VERDICT.txt" "$B/BGC_SUMMARY.tsv" "$B/bgc_regions.tsv"

# antiSMASH output directories are the expensive thing here (minutes per genome)
# and are the cache, so they are NOT wiped by default. What IS wiped is any
# directory that contains no usable result — failed-run debris, which would
# otherwise abort the retry with "Output directory contains other files".
NDEB=0
for gd in "$B"/as/*/; do
  [ -d "$gd" ] || continue
  gd="${gd%/}"
  if ls "$gd"/*.json >/dev/null 2>&1 || [ -s "$gd/index.html" ]; then continue; fi
  case "$gd" in "$B/as/"?*) rm -rf "${gd:?}"; NDEB=$((NDEB+1)) ;; esac
done
[ "$NDEB" -gt 0 ] && echo "  preflight: removed $NDEB empty/failed antiSMASH output dir(s)"

# FRESH=1 forces a complete rebuild: every antiSMASH result, every cached HMM
# search, every stamp. Use it when the METHOD changed, not to clear an error —
# a re-run already clears anything that is not a result.
if [ "${FRESH:-0}" = "1" ]; then
  echo "  FRESH=1 — discarding ALL cached results for this step"
  case "$B" in
    */33_antismash) rm -rf "${B:?}/as" "${B:?}/hmm" "${B:?}/as_logs" ;;
    *) echo "  !! refusing to clear '$B' — unexpected path"; exit 1 ;;
  esac
  mkdir -p "$B/as" "$B/hmm"
fi

BREADTH_INSITU="${BREADTH_INSITU:-0.50}"   # step 26's presence threshold
MEGA_LEN="${MEGA_LEN:-1000}"               # aa; megasynthase floor

echo "######## 33 — BIOSYNTHETIC GENE CLUSTERS ########"
echo "  started: $(date -Is)"
runlog "33_antismash started"

NGEN=$(ls -1d "$A"/*/ 2>/dev/null | wc -l)
echo "  genomes with a Bakta annotation: $NGEN"
if [ "$NGEN" -lt 1 ]; then echo "  !! no annotations in $A — step 10 has not run"; exit 1; fi

# =============================================================================
# 0. IS antiSMASH ACTUALLY USABLE?
#    Installed is not the same as usable. antiSMASH without its databases runs
#    and produces empty results, which would read as "no BGCs" — the worst
#    possible failure mode for this step. Both checks must pass.
# =============================================================================
echo
echo "######## 0. LOOKING FOR antiSMASH ########"

# WHY CONDA CANNOT INSTALL antiSMASH HERE, AND WHY THAT IS NOT A LOCAL PROBLEM
#   `conda create -n antismash antismash` fails on this machine, and the solver
#   output says why: every antiSMASH build pins an exact biopython (1.78, 1.81,
#   <=1.79 ...) and the pre-8 builds drag in ete2 -> mysql-python, which is a
#   Python 2 package. Those pins conflict with what is already resolvable in
#   this channel set. This is antiSMASH's dependency graph, not a broken env,
#   and no amount of retrying conda will change it. The project's own answer is
#   the container, which carries a resolved environment AND the databases.
#
# THREE TIERS, IN ORDER
#   A1  native antiSMASH in a conda env, if one exists and its databases check
#   A2  antiSMASH in Docker/Podman via the official wrapper  <-- the route here
#   B   the Pfam core-domain floor
AS_ENV=""; AS_VER=""; AS_OK=0; AS_HOW=""
# Same reasoning as resolve_hmm_tool below: `conda run` does not work on this
# machine, so a native antiSMASH is invoked through its env prefix when the
# patched config has given us env_prefix().
as_run() {                       # as_run <env> <args...>
  local e="$1"; shift
  local prefix=""
  if [ "$e" = "-" ]; then antismash "$@"; return $?; fi
  if declare -f env_prefix >/dev/null 2>&1; then
    prefix=$(env_prefix "$e" 2>/dev/null) || prefix=""
  fi
  if [ -n "$prefix" ] && [ -x "$prefix/bin/antismash" ]; then
    PATH="$prefix/bin:$PATH" LD_LIBRARY_PATH="$prefix/lib:${LD_LIBRARY_PATH:-}" \
      CONDA_PREFIX="$prefix" "$prefix/bin/antismash" "$@"
  else
    conda run --no-capture-output -n "$e" antismash "$@"
  fi
}

# ------------------------------------------------------------- A1: native conda
CANDS=()
[ -n "${ENV_antismash:-}" ] && CANDS+=("$ENV_antismash")
command -v antismash >/dev/null 2>&1 && CANDS+=("-")
while read -r e; do [ -n "$e" ] && CANDS+=("$e"); done < <(
  conda env list 2>/dev/null | awk '{print $1}' | grep -iE 'antismash|antiSMASH|secmet' || true)

for e in "${CANDS[@]:-}"; do
  [ -n "$e" ] || continue
  v=$(as_run "$e" --version 2>/dev/null | head -1)
  if [ -n "$v" ]; then AS_ENV="$e"; AS_VER="$v"; break; fi
done

if [ -n "$AS_ENV" ]; then
  echo "  native install found: $AS_VER   (env: $AS_ENV)"
  echo "  checking databases with --check-prereqs ..."
  if as_run "$AS_ENV" --check-prereqs > "$B/as/prereqs.txt" 2>&1; then
    AS_OK=1; AS_HOW="native"
    echo "  [PASS] prerequisites and databases satisfied"
    record_version "antismash" "$AS_ENV" antismash --version
  else
    echo "  [FAIL] --check-prereqs did not pass. Last lines:"
    tail -12 "$B/as/prereqs.txt" | sed 's/^/         /'
    echo "         Installed but its databases are missing. antiSMASH would"
    echo "         still RUN and return zero clusters, and a fabricated zero is"
    echo "         worse here than no antiSMASH at all. Trying the container."
  fi
else
  echo "  no native antiSMASH in any conda environment."
fi

# ------------------------------------------------------------ A2: the container
# The standalone-FULL image ships the databases inside it, so there is nothing
# to download separately and nothing to point at — which also means there is no
# way for this step to silently run against a half-populated database, the
# failure mode the prereq check above exists to catch.
CE=""; WRAP="${ANTISMASH_WRAPPER:-$WS2/bin/run_antismash}"
if [ "$AS_OK" != "1" ]; then
  echo
  echo "  --- container route ---"
  for c in docker podman; do
    command -v "$c" >/dev/null 2>&1 || continue
    if "$c" info >/dev/null 2>&1; then CE="$c"; break; fi
    echo "  $c is installed but '$c info' failed — the daemon is not reachable"
    echo "  as this user. For docker that is usually group membership:"
    echo "      sudo usermod -aG docker \$USER   # then log out and back in"
  done

  if [ -z "$CE" ]; then
    echo "  no usable container engine (tried docker, podman)."
  else
    echo "  container engine: $CE ($("$CE" --version 2>/dev/null | head -1))"

    # Disk. The full image is ~9 GB and the engine stores it under its own root
    # (/var/lib/docker for docker), which is on '/' here — the filesystem this
    # project is under standing orders to keep clear. Check before pulling
    # rather than discovering it half way through a 9 GB download.
    CROOT=$("$CE" info --format '{{.DockerRootDir}}' 2>/dev/null)
    [ -n "$CROOT" ] || CROOT="/var/lib/${CE}"
    AVAILG=$(df -BG --output=avail "$CROOT" 2>/dev/null | tail -1 | tr -dc '0-9')
    AVAILG=${AVAILG:-0}
    echo "  image store: $CROOT   free: ${AVAILG} GB   (image needs ~9 GB)"
    if [ "$AVAILG" -lt 15 ]; then
      echo "  !! under 15 GB free where $CE keeps images. Not pulling."
      echo "     Move the image store to ${SCRATCH_ROOT} first — for docker,"
      echo "     set  {\"data-root\": \"${SCRATCH_ROOT}/docker\"}  in"
      echo "     /etc/docker/daemon.json and restart the daemon."
      CE=""
    fi
  fi

  if [ -n "$CE" ]; then
    mkdir -p "$WS2/bin"
    if [ ! -x "$WRAP" ]; then
      echo "  fetching the official wrapper script ..."
      curl -fsSL -o "$WRAP.part" \
        https://dl.secondarymetabolites.org/releases/latest/docker-run_antismash-full \
        2>"$B/as/curl.err"
      # A wrapper that is actually an HTML error page would run and fail in a
      # confusing way later. Same lesson as the dbCAN download at step 30.
      if [ -s "$WRAP.part" ] && head -c 2 "$WRAP.part" | grep -q '#!'; then
        mv "$WRAP.part" "$WRAP"; chmod +x "$WRAP"
        echo "  wrapper: $WRAP"
      else
        rm -f "$WRAP.part"
        echo "  !! download did not return a shell script."
        [ -s "$B/as/curl.err" ] && sed 's/^/     /' "$B/as/curl.err"
        CE=""
      fi
    else
      echo "  wrapper already present: $WRAP"
    fi
  fi

  # The wrapper is written for docker. The upstream instructions say to swap the
  # word for podman; do that on a COPY so the original is never edited in place.
  if [ -n "$CE" ] && [ "$CE" = "podman" ]; then
    sed 's/\bdocker\b/podman/g' "$WRAP" > "$WS2/bin/run_antismash.podman"
    chmod +x "$WS2/bin/run_antismash.podman"; WRAP="$WS2/bin/run_antismash.podman"
    echo "  podman variant: $WRAP"
  fi

  # The image name is READ OUT OF THE WRAPPER, not written here from memory.
  # The wrapper is the upstream artefact that knows which image goes with this
  # release; hardcoding a tag would be exactly the kind of asserted fact this
  # project keeps getting wrong.
  IMG=""; CDST=""
  if [ -n "$CE" ]; then
    IMG=$(grep -oE '(docker\.io/)?antismash/[a-zA-Z0-9._-]+' "$WRAP" 2>/dev/null | head -1)
    # The wrapper pins a version in a variable and interpolates it into the
    # image reference, so the tag has to be read separately or we end up pulling
    # :latest while the wrapper runs :8.0.4 — two different programs.
    ASV=$(sed -n 's/^readonly VERSION="\([^"]*\)".*/\1/p' "$WRAP" 2>/dev/null | head -1)
    [ -n "$IMG" ] && [ -n "$ASV" ] && IMG="${IMG}:${ASV}"
    # Where the wrapper mounts the output directory INSIDE the container. Read
    # it rather than assume it: antiSMASH needs to be told to write there.
    CDST=$(sed -n 's/^readonly CONTAINER_DST_DIR=//p' "$WRAP" 2>/dev/null | head -1)
    CDST="${CDST%\"}"; CDST="${CDST#\"}"
    if [ -z "$IMG" ]; then
      echo "  !! could not read an image name out of the wrapper. Inspect it:"
      echo "     less $WRAP"
      CE=""
    elif [ -z "$CDST" ]; then
      echo "  !! could not read CONTAINER_DST_DIR out of the wrapper."
      echo "     Without it antiSMASH cannot be told where to write. Inspect:"
      echo "     head -40 $WRAP"
      CE=""
    else
      echo "  image (read from the wrapper): $IMG"
      echo "  container output path: $CDST"
      echo "  pulling — the full image is ~9 GB, so the first run is slow ..."
      if "$CE" pull "$IMG" 2>&1 | tail -3 | sed 's/^/     /'; then :; fi
      if ! "$CE" image inspect "$IMG" >/dev/null 2>&1; then
        echo "  !! image not present after pull. Nothing to run."
        CE=""
      fi
    fi
  fi

  if [ -n "$CE" ]; then
    # Read the version out of the container. Two spellings are tried because
    # whether the image's entrypoint IS antismash or a shell differs between
    # releases; whichever answers is the truth for this image.
    : > "$B/as/version.txt"
    "$CE" run --rm "$IMG" --version           >>"$B/as/version.txt" 2>&1 || true
    "$CE" run --rm "$IMG" antismash --version >>"$B/as/version.txt" 2>&1 || true
    AS_VER=$(grep -oiE 'antismash[ -]?[0-9][0-9.]*' "$B/as/version.txt" | head -1)
    AS_OK=1; AS_HOW="container"
    if [ -n "$AS_VER" ]; then
      echo "  [PASS] $AS_VER (via $CE)"
    else
      AS_VER="unknown (image $IMG)"
      echo "  [WARN] the image is present but would not print a version."
      echo "         Proceeding — the per-genome runs below will show whether it"
      echo "         works — but the version table records it as unknown, and a"
      echo "         thesis methods section needs the real number. Get it with:"
      echo "           $CE run --rm $IMG --version"
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' "$STEP" "antismash" "$AS_VER" \
      "$CE $IMG" "$(date -Is)" >> "$VERSIONS"
  fi
fi

# --------------------------------------------------------- CAN IT SEE OUR DISK?
# A container that cannot bind-mount the workspace still RUNS and still EXITS 0.
# It reads nothing, writes into its own throwaway layer, and reports success.
# That is what happened on the last attempt: 29 genomes "completed" in seconds
# and not one output file reached the disk.
#
# Docker here is the SNAP package (image store under /var/snap/docker/...), and
# snap confinement restricts which host paths the daemon may bind-mount. The
# workspace is on ${SCRATCH_ROOT}, which is outside anything snap grants by
# default. So the mount silently yields nothing.
#
# This probe settles it before 29 genomes are wasted: write a file on the host,
# read it from inside the container, write a file back out, and check it landed.
STAGE=""
probe_mount() {                  # probe_mount <hostdir>  -> 0 if usable
  local hd="$1"
  local p="$hd/.as_mount_probe.$$"
  mkdir -p "$p" 2>/dev/null || return 1
  echo "PROBE-IN" > "$p/in.txt" 2>/dev/null || { rm -rf "$p"; return 1; }
  # --user MUST match what the wrapper uses. Probing as root would pass on a
  # mount the real run cannot write to, which is exactly the failure we are
  # chasing: PermissionError on mkdir inside /output.
  # mkdir is included because that is the operation antiSMASH actually performs.
  "$CE" run --rm --entrypoint /bin/sh --user="$(id -u):$(id -g)" \
      -v "$p:/probe" "$IMG" \
      -c 'mkdir -p /probe/sub && cat /probe/in.txt > /probe/sub/out.txt' >/dev/null 2>&1
  local ok=1
  [ "$(cat "$p/sub/out.txt" 2>/dev/null)" = "PROBE-IN" ] && ok=0
  rm -rf "$p"
  return "$ok"
}

if [ "$AS_OK" = "1" ] && [ "$AS_HOW" = "container" ]; then
  echo
  echo "  --- can the container reach the workspace? ---"
  if probe_mount "$B"; then
    echo "  [PASS] bind mount works on the workspace filesystem."
  else
    echo "  [FAIL] the container cannot create a directory under $B."
    echo "         Symptom in the antiSMASH log:"
    echo "           PermissionError: [Errno 13] Permission denied: '$CDST/<genome>'"
    echo "         The bind mount is not delivering the host directory, so what"
    echo "         the container sees at $CDST is a root-owned empty directory"
    echo "         and --user=\$(id -u) cannot write into it."
    echo "         Looking for a directory it CAN reach ..."
    for cand in "$(readlink -f "$HOME" 2>/dev/null)" /tmp; do
      [ -n "$cand" ] && [ -d "$cand" ] || continue
      if probe_mount "$cand"; then
        STAGE="$cand/.antismash_stage"
        echo "  [PASS] staging through $cand"
        break
      fi
      echo "         $cand — no"
    done
    if [ -z "$STAGE" ]; then
      echo
      echo "  !! no directory on this machine is mountable by this Docker."
      echo "     Nothing is written. 29 zeros from a container that never saw the"
      echo "     data would be indistinguishable from 29 real negatives."
      echo
      echo "     This Docker is the SNAP package (image store under /var/snap),"
      echo "     and snap confinement decides which host paths its daemon may"
      echo "     bind-mount. Check what it currently has:"
      echo "       snap connections docker"
      echo
      echo "     Fixes, in order of least disruption:"
      echo "       1. sudo snap connect docker:removable-media"
      echo "          The removable-media interface covers /mnt, which is where"
      echo "          this workspace lives. If it works, nothing else changes and"
      echo "          no staging is needed. Try this first, then re-run."
      echo "       2. sudo snap connect docker:home"
      echo "          Grants \$HOME. The workspace is still on /mnt, so this step"
      echo "          would then stage through \$HOME automatically."
      echo "       3. Replace snap Docker with the apt package, which has no such"
      echo "          confinement:  sudo snap remove docker, then docker-ce from"
      echo "          docs.docker.com — set its data-root to ${SCRATCH_ROOT} so"
      echo "          the images do not fill /."
      echo "       4. podman, which is rootless and unconfined here."
      runlog "**33_antismash ABORTED: container cannot bind-mount host paths**"
      exit 6
    fi
    echo "         Genomes will be copied in and results copied back out."
    echo "         Cost: a few hundred MB of transient files under $STAGE."
    mkdir -p "$STAGE"
  fi
fi

if [ "$AS_OK" != "1" ]; then
  echo
  echo "  No usable antiSMASH by either route."
  echo "  Falling back to the Pfam core-domain floor (PATH B)."
fi

# antiSMASH runs through whichever route succeeded. Everything downstream is
# identical, so the route is recorded and then stops mattering.
# THE WRAPPER DOES NOT PASS --output-dir. It mounts the host output directory at
# CONTAINER_DST_DIR and then hands antiSMASH only the input filename. antiSMASH
# therefore falls back to its default, which is a subdirectory NAMED AFTER THE
# INPUT, inside the container's working directory:
#     PermissionError: [Errno 13] Permission denied: '/output/SED_das_semibin_002'
# Two consequences, and I had missed both. Results would land one level deeper
# than this script looks for them, and the mkdir needs write permission on the
# mount. Passing --output-dir with the CONTAINER-side path puts the files
# directly in the mounted directory, which is the layout the parser expects.
as_exec() {                      # as_exec <gbff> <outdir> <args...>
  local gb="$1"; shift
  local od="$1"; shift
  if [ "$AS_HOW" = "container" ]; then "$WRAP" "$gb" "$od" --output-dir "$CDST" "$@"
  else as_run "$AS_ENV" "$@" --output-dir "$od" "$gb"; fi
}

# =============================================================================
# 1. PATH A — antiSMASH
# =============================================================================
MODE="pfam"
if [ "$AS_OK" = "1" ]; then
MODE="antismash"
echo
echo "######## 1. RUNNING antiSMASH ########"
echo "  NOTE on --taxon: antiSMASH's rule sets are bacteria / fungi / plants."
echo "  There is no archaeal rule set, so the eight archaeal genomes are run"
echo "  under the bacterial rules. That is the only option available and it is"
echo "  stated rather than hidden: archaeal calls here are approximate."
echo
# antiSMASH refuses to write into a directory that already contains anything:
#   "Output directory contains other files, aborting for safety"
# The first version of this script created antismash.log INSIDE the output
# directory, so the redirection made the directory non-empty before antiSMASH
# had even started, and it aborted on all 29 genomes. Logs therefore live
# outside the output tree, and each output directory is emptied before use so a
# previous failure's debris cannot abort the retry either.
LOGD="$B/as_logs"; mkdir -p "$LOGD"

n=0; NAS_OK=0
for d in "$A"/*/; do
  name=$(basename "$d")
  gbff="$d/${name}.gbff"
  [ -s "$gbff" ] || { echo "  [skip] $name — no gbff"; continue; }
  n=$((n+1))
  out="$B/as/$name"
  # A genome with genuinely zero clusters still produces a JSON and an
  # index.html; only the region GBKs are absent. Keying the cache on region
  # files alone would re-run every clean negative forever.
  if ls "$out"/*.json >/dev/null 2>&1 || [ -s "$out/index.html" ]; then
    printf '  [%2d/%2d] %-30s cached\n' "$n" "$NGEN" "$name"
    NAS_OK=$((NAS_OK+1)); continue
  fi
  printf '  [%2d/%2d] %-30s running ...' "$n" "$NGEN" "$name"
  log="$LOGD/${name}.log"

  # Where antiSMASH is actually pointed. With staging this is under $STAGE and
  # the results are copied back afterwards; without it, it is the final home.
  if [ -n "$STAGE" ]; then
    wdir="$STAGE/$name"; win="$STAGE/in_${name}.gbff"
    rm -rf "$wdir"; mkdir -p "$wdir"; cp -f "$gbff" "$win"
  else
    wdir="$out"; win="$gbff"
  fi

  # Empty the output directory. Guarded on the path shape so a mistake in $B or
  # $name can never turn this into a recursive delete somewhere else.
  clear_out() {
    case "$out" in
      "$B/as/"?*) rm -rf "${out:?}" ;;
      *) echo " !! refusing to clear '$out' — unexpected path"; return 1 ;;
    esac
    mkdir -p "$out"
    [ -n "$STAGE" ] && { rm -rf "$wdir"; mkdir -p "$wdir"; }
    return 0
  }
  clear_out || { echo; continue; }

  as_exec "$win" "$wdir" \
      --genefinding-tool none \
      --taxon bacteria \
      --cb-general --cb-knownclusters --cb-subclusters \
      --asf --pfam2go \
      --cpus "$THREADS" \
      --output-basename "$name" > "$log" 2>&1
  rc=$?
  # An option this build does not know is a hard exit before any analysis runs.
  # Retry once with the minimum that every version from 5 onward accepts, and
  # SAY that the extra analyses were dropped — a quieter run is not the same run.
  if [ "$rc" -ne 0 ] && grep -qiE 'unrecognized arguments|no such option|invalid choice' \
        "$log" 2>/dev/null; then
    echo -n " retrying without optional modules ..."
    touch "$LOGD/${name}.REDUCED"
    clear_out || { echo; continue; }
    as_exec "$win" "$wdir" --genefinding-tool none --taxon bacteria \
        --cpus "$THREADS" >> "$log" 2>&1
    rc=$?
  fi

  if [ -n "$STAGE" ]; then
    cp -a "$wdir/." "$out/" 2>/dev/null
    rm -rf "$wdir" "$win"
  fi

  # SUCCESS IS OUTPUT ON DISK, NOT EXIT ZERO.
  # The previous run took rc=0 as completion. All 29 exited 0 having produced
  # nothing, NAS_OK reached 29, the "did anything run?" guard was satisfied, and
  # a table of 29 zeros passed every validation check. An exit code describes
  # the process; only the files describe the analysis.
  if ls "$out"/*.json >/dev/null 2>&1 || [ -s "$out/index.html" ]; then
    nreg=$(ls -1 "$out"/*.region*.gbk 2>/dev/null | wc -l)
    echo " done — $nreg region(s)"
    NAS_OK=$((NAS_OK+1))
  elif [ "$rc" -eq 0 ]; then
    echo " EXITED 0 BUT WROTE NOTHING"
    echo "           No JSON and no index.html in $out — the run did not happen."
    tail -6 "$log" | sed 's/^/           /'
  else
    echo " FAILED (rc=$rc)"
    tail -8 "$log" | sed 's/^/           /'
  fi
done
[ -n "$STAGE" ] && rm -rf "$STAGE"

# A genome that never ran contributes a zero to every count below, and zeros
# from a failed run are indistinguishable from zeros from a clean negative once
# they are in a table. That is the shape of the step-25 provirus bug. Refuse.
echo
echo "  genomes antiSMASH completed: $NAS_OK/$n"
if [ "$NAS_OK" -eq 0 ]; then
  echo "  !! not one genome completed. Nothing is written."
  echo "     'No BGCs in any genome' and 'antiSMASH never ran' produce the same"
  echo "     table, and only one of them is a result."
  runlog "**33_antismash ABORTED: 0/$n genomes completed**"
  exit 5
elif [ "$NAS_OK" -lt "$n" ]; then
  echo "  !! $((n - NAS_OK)) genome(s) failed. They appear below with zero"
  echo "     regions, which for them means NOT MEASURED, not absent."
fi
NRED=$(ls -1 "$LOGD"/*.REDUCED 2>/dev/null | wc -l)
if [ "$NRED" -gt 0 ]; then
  echo
  echo "  !! $NRED genome(s) ran with the reduced option set: no ClusterBlast,"
  echo "     no MIBiG comparison, no active-site finder. Cluster CALLS are"
  echo "     unaffected — those come from the core rules — but similarity to"
  echo "     known clusters is missing for those genomes and must not be"
  echo "     reported as 'no known cluster'."
fi

# ------------------------------------------------- parse regions out of the run
echo
echo "  parsing regions ..."
python3 - "$B" <<'PYA'
import sys, os, glob, json, csv, re

B = sys.argv[1]
rows = []
for d in sorted(glob.glob(os.path.join(B, 'as', '*'))):
    if not os.path.isdir(d): continue
    name = os.path.basename(d)

    # Preferred: the JSON summary. Fall back to the per-region GenBank files,
    # whose 'region' feature carries /product= . Both are antiSMASH's own
    # output; neither is a guess about what it meant.
    got = False
    js = os.path.join(d, f'{name}.json')
    if not os.path.exists(js):
        # the reduced-option retry drops --output-basename, so the JSON is
        # named after the input record instead
        cands = [p for p in glob.glob(os.path.join(d, '*.json'))
                 if not p.endswith('regions.js')]
        js = cands[0] if cands else js
    if os.path.exists(js):
        try:
            data = json.load(open(js))
            for rec in data.get('records', []):
                cid = rec.get('id', '?')
                for ar in rec.get('areas', []):
                    prods = ar.get('products') or []
                    rows.append(dict(genome=name, contig=cid,
                                     start=ar.get('start', ''), end=ar.get('end', ''),
                                     products=';'.join(prods) or 'unknown',
                                     source='json'))
                    got = True
        except Exception as e:
            print(f"    [warn] {name}: JSON unreadable ({e}) — using region GBKs")

    if not got:
        # CAREFUL: /product= appears on EVERY CDS in a region GenBank as well as
        # on the region feature itself. Collecting all of them would fill the
        # product column with "hypothetical protein" and call it a cluster type.
        # GenBank puts a feature key at column 6, so a line matching five spaces
        # then a non-space starts a new feature; only /product lines while that
        # feature is 'region' (or 'protocluster') are cluster products.
        FEATKEY = re.compile(r'^ {5}(\S+)\s')
        for g in sorted(glob.glob(os.path.join(d, '*.region*.gbk'))):
            prods, cid, inreg = [], os.path.basename(g), False
            for line in open(g, errors='replace'):
                if line.startswith('LOCUS'):
                    parts = line.split()
                    if len(parts) > 1: cid = parts[1]
                    continue
                m = FEATKEY.match(line)
                if m:
                    inreg = m.group(1) in ('region', 'protocluster', 'cand_cluster')
                    continue
                s = line.strip()
                if inreg and s.startswith('/product='):
                    prods.append(s.split('=', 1)[1].strip('"'))
            rows.append(dict(genome=name, contig=cid, start='', end='',
                             products=';'.join(sorted(set(prods))) or 'unknown',
                             source='gbk'))

with open(os.path.join(B, 'bgc_regions.tsv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, delimiter='\t',
                       fieldnames=['genome', 'contig', 'start', 'end', 'products', 'source'])
    w.writeheader(); w.writerows(rows)
print(f"    {len(rows)} region(s) across {len({r['genome'] for r in rows})} genome(s)")
PYA
fi

# =============================================================================
# 2. PATH B — Pfam core-domain floor
# =============================================================================
if [ "$MODE" = "pfam" ]; then
echo
echo "######## 1. PFAM CORE-DOMAIN FLOOR ########"

validate_hmm() {   # a downloaded HTML error page is not a database (step 30)
  local f="$1"
  [ -s "$f" ] || return 1
  head -c 8 "$f" 2>/dev/null | grep -q '^HMMER3' || return 1
  return 0
}

HMMDB=""; HMMBASE=""; HMODE=""
for c in "${PFAM_HMM:-}" "$BAKTA_DB/pfam.hmm" "$BAKTA_DB/Pfam-A.hmm" \
         "$WS2/00_db/pfam/Pfam-A.hmm" \
         ${DB_DIR}/pfam/Pfam-A.hmm; do
  [ -n "$c" ] && validate_hmm "$c" && { HMMDB="$c"; HMODE="search"; break; }
done
if [ -z "$HMODE" ]; then
  for h3 in "${PFAM_H3M:-}" "$BAKTA_DB"/pfam.h3m "$BAKTA_DB"/Pfam-A.h3m \
            $(find "$BAKTA_DB" ${SCRATCH_ROOT} -maxdepth 5 \
                   \( -name 'pfam*.h3m' -o -name 'Pfam*.h3m' \) 2>/dev/null | head -5); do
    [ -n "${h3:-}" ] && [ -s "$h3" ] || continue
    base="${h3%.h3m}"
    ok=1; for ext in h3i h3f h3p; do [ -s "${base}.${ext}" ] || ok=0; done
    [ "$ok" = "1" ] && { HMMBASE="$base"; HMODE="scan"; break; }
  done
fi

if [ -z "$HMODE" ]; then
  echo
  echo "  !! No usable Pfam HMM library and no working antiSMASH."
  echo
  echo "     Nothing is written. An empty BGC table would be read as 'this"
  echo "     community makes no secondary metabolites', which is a claim this"
  echo "     run has not tested."
  echo
  echo "     Either of these fixes it:"
  echo "       1. antiSMASH, which is the analysis this step actually wants:"
  echo "            conda create -n antismash -c bioconda -c conda-forge antismash"
  echo "            conda run -n antismash download-antismash-databases"
  echo "          (~10 GB; put it on ${SCRATCH_ROOT}, not /)"
  echo "       2. Point at a Pfam library already on disk:"
  echo "            PFAM_HMM=/path/to/Pfam-A.hmm  bash scripts/33_antismash_bgc.sh"
  echo "            PFAM_H3M=/path/to/pfam.h3m    bash scripts/33_antismash_bgc.sh"
  echo "          The Bakta database normally ships one — checked $BAKTA_DB."
  runlog "33_antismash ABORTED: no antiSMASH and no Pfam library"
  exit 3
fi

# ------------------------------------------------- select models by their DESC
# In 'search' mode the library is plain text, so DESC is read here and used to
# build a small subset — fast. In 'scan' mode the library is binary, so the whole
# library is scanned and DESC is recovered from hmmscan's own output. Same
# selection criterion either way, and in neither case is an accession asserted.
SUBSET="$B/hmm/core_subset.hmm"
NAMEMAP="$B/hmm/model_desc.tsv"

if [ "$HMODE" = "search" ]; then
  echo "  library: $HMMDB (plain text)"
  DBHASH=$(md5sum "$HMMDB" | cut -c1-12)
  python3 - "$HMMDB" "$NAMEMAP" <<'PYS'
import sys, re
src, out = sys.argv[1], sys.argv[2]
name = acc = desc = None
rows = []
for line in open(src, errors='replace'):
    if line.startswith('NAME'):  name = line.split(None, 1)[1].strip(); acc = desc = None
    elif line.startswith('ACC'): acc  = line.split(None, 1)[1].strip()
    elif line.startswith('DESC'):desc = line.split(None, 1)[1].strip()
    elif line.startswith('//'):
        if name: rows.append((name, acc or '-', desc or ''))
        name = acc = desc = None
with open(out, 'w') as fh:
    fh.write('model\taccession\tdesc\n')
    for r in rows: fh.write('\t'.join(r) + '\n')
print(f"  models in library: {len(rows)}")
PYS
  NMOD=$(($(wc -l < "$NAMEMAP") - 1))
  if [ "$NMOD" -lt 1000 ]; then
    echo "  !! only $NMOD models parsed — that does not look like Pfam"; exit 1
  fi
else
  echo "  library: ${HMMBASE}.h3m (pressed — DESC comes from hmmscan output)"
  DBHASH=$(md5sum "${HMMBASE}.h3m" | cut -c1-12)
  : > "$NAMEMAP"
fi
echo "  library hash: $DBHASH"

# The keyword patterns. These are SEARCH TERMS against Pfam DESC lines, not
# assertions about which families exist. Every model each pattern actually
# matched is printed below.
cat > "$B/hmm/patterns.tsv" <<'PAT'
class	confidence	pattern
NRPS_C	high	condensation domain|condensation protein
NRPS_A	shared	amp-binding|aminoacid adenylation|amino acid adenylation|adenylation domain
CARRIER	shared	phosphopantetheine|acyl carrier protein
TE	shared	thioesterase
PKS_KS	shared	beta-ketoacyl synth|ketoacyl synthase|ketoacyl-synt|ketosynthase
PKS_AT	shared	acyl transferase domain|acyltransferase domain|malonyl.*transacylase
PKS_KR	shared	ketoreductase|ketoacyl reductase
TERPENE	medium	squalene.*synth|phytoene synth|terpene synth|terpene cyclase|prenyltransferase.*aromatic
RIPP_LANTI	high	lanthionine|lantibiotic|lanthipeptide
RIPP_YCAO	medium	ycao|cyclodehydratase|thiazole.*oxazole
RIPP_SPASM	medium	spasm domain|twitch domain
BACTERIOCIN	high	bacteriocin|microcin|colicin|lasso peptide|sactipeptide
SIDEROPHORE	high	siderophore|iuca|iucc|nis synthetase|aerobactin
PHOSPHONATE	high	phosphoenolpyruvate.*mutase|phosphonopyruvate|phosphonate biosynth
BETALACTONE	medium	beta-lactone|betalactone
PAT

if [ "$HMODE" = "search" ]; then
  python3 - "$NAMEMAP" "$B/hmm/patterns.tsv" "$B/hmm/selected.tsv" <<'PYX'
import sys, csv, re
nm, pf, out = sys.argv[1], sys.argv[2], sys.argv[3]
pats = [(r['class'], r['confidence'], re.compile(r['pattern'], re.I))
        for r in csv.DictReader(open(pf), delimiter='\t')]
sel = []
for r in csv.DictReader(open(nm), delimiter='\t'):
    d = r['desc'] or ''
    n = r['model'] or ''
    for cls, conf, rx in pats:
        if rx.search(d) or rx.search(n):
            sel.append((r['model'], r['accession'], cls, conf, d)); break
with open(out, 'w') as fh:
    fh.write('model\taccession\tclass\tconfidence\tdesc\n')
    for s in sel: fh.write('\t'.join(s) + '\n')
from collections import Counter
c = Counter(s[2] for s in sel)
print(f"  selected {len(sel)} models:")
for k in sorted(c): print(f"    {k:<14} {c[k]:>4}")
PYX
  NSEL=$(($(wc -l < "$B/hmm/selected.tsv") - 1))
  if [ "$NSEL" -lt 10 ]; then
    echo "  !! only $NSEL models matched the patterns — the DESC parse or the"
    echo "     patterns are wrong. Inspect $NAMEMAP before trusting anything."
    exit 1
  fi
  echo
  echo "  --- every selected model, for audit ---"
  column -t -s$'\t' "$B/hmm/selected.tsv" 2>/dev/null | cut -c1-150 | head -80
  echo

  # Build the subset library by copying whole model blocks out of the plain-text
  # file. NOT hmmfetch: hmmfetch wants an .ssi index written NEXT TO the library,
  # and the Bakta database directory is not guaranteed writable — the same class
  # of assumption that broke NB_Genome's conda metadata. A model block is
  # everything from its HMMER3 header line to the '//' terminator, so copying it
  # verbatim preserves the GA thresholds, which is what --cut_ga needs.
  cut -f1 "$B/hmm/selected.tsv" | tail -n +2 > "$B/hmm/names.txt"
  python3 - "$HMMDB" "$B/hmm/names.txt" "$SUBSET" <<'PYF'
import sys
src, names, out = sys.argv[1], sys.argv[2], sys.argv[3]
want = {l.strip() for l in open(names) if l.strip()}
block, cur, kept = [], None, 0
with open(out, 'w') as fh:
    for line in open(src, errors='replace'):
        block.append(line)
        if line.startswith('NAME'):
            cur = line.split(None, 1)[1].strip()
        elif line.startswith('//'):
            if cur in want:
                fh.writelines(block); kept += 1
            block, cur = [], None
print(f"  subset library: {kept}/{len(want)} model(s) copied")
if kept != len(want):
    print(f"  !! {len(want)-kept} selected model(s) were not found in the library")
PYF
  if ! validate_hmm "$SUBSET"; then
    echo "  !! subset library is not a HMMER file — aborting rather than searching"
    echo "     against something unreadable and reporting the resulting zero."
    exit 1
  fi
fi

# --------------------------------------------------- resolve the HMMER binary
# The first run of this step died 29 times with "no working env resolved for
# hmmscan". env_resolved.sh carries ENV_hmmsearch because that is the binary
# the smoketest exercised — but hmmsearch and hmmscan ship in the SAME hmmer
# package, so the env that has one has the other. run_tool cannot know that: it
# looks up ENV_<tool> and gives up.
#
# The fix is not to assume it. Each candidate env is tested by actually running
# the binary, and only an env that answers is used. An assumed env that turns
# out to be wrong would produce an empty domtbl, and an empty domtbl is a zero.
#
# NOTE ON HOW THE TEST IS DONE. 00_config_patch.sh replaced run_tool so it stops
# using `conda run` entirely — conda's metadata for NB_Genome is damaged and the
# directory is not writable, so `conda run -n NB_Genome` fails even though the
# binaries are fine. Testing candidates with `conda run` would therefore reject
# every environment on this machine. The test below uses env_prefix() when the
# patch has defined it, and only falls back to conda run on an unpatched config.
resolve_hmm_tool() {             # resolve_hmm_tool <tool>  -> prints env or ''
  local tool="$1"
  local var="ENV_${tool}"
  local cands="${!var:-} ${ENV_hmmsearch:-} ${ENV_hmmscan:-} ${ENV_hmmer:-} ${ENV_bakta:-} -"
  local e prefix
  for e in $cands; do
    [ -n "$e" ] || continue
    if [ "$e" = "-" ]; then
      command -v "$tool" >/dev/null 2>&1 && { echo "-"; return 0; }
      continue
    fi
    if declare -f env_prefix >/dev/null 2>&1; then
      prefix=$(env_prefix "$e" 2>/dev/null) || continue
      [ -x "$prefix/bin/$tool" ] && { echo "$e"; return 0; }
    else
      conda run -n "$e" "$tool" -h >/dev/null 2>&1 && { echo "$e"; return 0; }
    fi
  done
  return 1
}

for t in hmmsearch hmmscan; do
  v="ENV_${t}"
  cur="${!v:-}"
  got=$(resolve_hmm_tool "$t" || true)
  if [ -z "$got" ]; then
    echo "  !! $t does not run in any candidate environment."
    echo "     Tried: ENV_$t='${cur}', ENV_hmmsearch='${ENV_hmmsearch:-}',"
    echo "     ENV_hmmer='${ENV_hmmer:-}', ENV_bakta='${ENV_bakta:-}', and PATH."
    echo "     Bakta needs HMMER, so the env that runs bakta should have it."
    echo "     Nothing is written — a search that never ran is not a zero."
    runlog "33_antismash ABORTED: $t unresolvable"
    exit 4
  fi
  export "$v=$got"
  echo "  $t -> env '$got'$([ "$got" != "$cur" ] && echo '  (resolved here, not in env_resolved.sh)')"
done

# ------------------------------------------------------------------ the search
echo
echo "######## 2. SEARCHING $NGEN PROTEOMES ########"
if [ "$HMODE" = "scan" ]; then
  echo "  mode: hmmscan over the FULL library (pressed-only, and DESC is needed)."
  echo "  This is the slow direction. Expect 1-3 min per genome."
else
  echo "  mode: hmmsearch over the selected subset."
fi

n=0; nfail=0
for d in "$A"/*/; do
  name=$(basename "$d"); faa="$d/${name}.faa"
  [ -s "$faa" ] || continue
  n=$((n+1))
  out="$B/hmm/${name}.domtbl"; stamp="$B/hmm/${name}.stamp"
  if [ -s "$out" ] && [ "$(cat "$stamp" 2>/dev/null)" = "${DBHASH}_${HMODE}" ]; then
    printf '\r  [%2d/%2d] %-30s cached   ' "$n" "$NGEN" "$name"; continue
  fi
  printf '\r  [%2d/%2d] %-30s searching' "$n" "$NGEN" "$name"
  if [ "$HMODE" = "search" ]; then
    run_tool hmmsearch --cut_ga --domtblout "$out.part" --cpu "$THREADS" \
        "$SUBSET" "$faa" >/dev/null 2>"$B/hmm/${name}.err"
    rc=$?
    if [ "$rc" -ne 0 ]; then     # library may carry no GA lines
      run_tool hmmsearch -E 1e-5 --domE 1e-5 --domtblout "$out.part" --cpu "$THREADS" \
          "$SUBSET" "$faa" >/dev/null 2>"$B/hmm/${name}.err"
      rc=$?
    fi
  else
    run_tool hmmscan --cut_ga --domtblout "$out.part" --cpu "$THREADS" \
        "$HMMBASE" "$faa" >/dev/null 2>"$B/hmm/${name}.err"
    rc=$?
    if [ "$rc" -ne 0 ]; then
      run_tool hmmscan -E 1e-5 --domE 1e-5 --domtblout "$out.part" --cpu "$THREADS" \
          "$HMMBASE" "$faa" >/dev/null 2>"$B/hmm/${name}.err"
      rc=$?
    fi
  fi
  if [ "$rc" -eq 0 ] && [ -s "$out.part" ]; then
    mv "$out.part" "$out"; echo "${DBHASH}_${HMODE}" > "$stamp"
  else
    nfail=$((nfail+1)); rm -f "$out.part"
    printf '\r  [%2d/%2d] %-30s FAILED   \n' "$n" "$NGEN" "$name"
    tail -3 "$B/hmm/${name}.err" 2>/dev/null | sed 's/^/           /'
  fi
done
printf '\r  %-64s\n' "  searched $n proteome(s), $nfail failure(s)"
if [ "$nfail" -gt 0 ]; then
  echo "  !! $nfail genome(s) have no result. Any zero for those is a missing"
  echo "     measurement, not an absence, and they are excluded below."
fi
fi   # end PATH B search

# =============================================================================
# 3. INTERPRETATION
# =============================================================================
echo
python3 - "$B" "$A" "$F" "$X" "$MODE" "${HMODE:-none}" "$BREADTH_INSITU" "$MEGA_LEN" <<'PY'
import sys, os, csv, glob, re
from collections import defaultdict, Counter

B, A, F, X, MODE, HMODE, BINSITU, MEGA = sys.argv[1:9]
BINSITU = float(BINSITU); MEGA = int(MEGA)
W = 108

def rule(t=''):
    print('=' * W)
    if t: print('  ' + t); print('=' * W)

# ------------------------------------------------------------------- taxonomy
tax = {}
tpath = os.path.join(F, 'MAG_TAXONOMY.tsv')
if os.path.exists(tpath):
    for r in csv.DictReader(open(tpath), delimiter='\t'):
        wn = r.get('working_name') or r.get('mag_id')
        if wn: tax[wn] = r
def g(m, k, d=''): return (tax.get(m, {}) or {}).get(k, '') or d
def mid(m):   return g(m, 'mag_id', m)
def genus(m): return g(m, 'genus', '-')
def phylum(m):return g(m, 'phylum', '')

# Domain: use an explicit column if there is one, then the GTDB classification
# string, and only then a phylum list. Which route was used is printed, because
# the bacterial control below is only valid if archaea were correctly excluded.
ARCH_PHYLA = {'halobacteriota', 'methanobacteriota', 'thermoplasmatota',
              'euryarchaeota', 'crenarchaeota', 'thermoproteota',
              'asgardarchaeota', 'nanoarchaeota', 'micrarchaeota'}
def domain_of(m):
    r = tax.get(m, {}) or {}
    for k in ('domain', 'Domain'):
        v = (r.get(k) or '').strip()
        if v: return 'Archaea' if 'archae' in v.lower() else 'Bacteria'
    cl = (r.get('classification') or r.get('gtdb_classification') or '')
    if 'd__Archaea' in cl:  return 'Archaea'
    if 'd__Bacteria' in cl: return 'Bacteria'
    return 'Archaea' if phylum(m).strip().lower() in ARCH_PHYLA else 'Bacteria'

route = 'phylum list (fallback)'
if tax:
    s = next(iter(tax.values()))
    if any(k in s for k in ('domain', 'Domain')): route = 'explicit domain column'
    elif any(k in s for k in ('classification', 'gtdb_classification')): route = 'GTDB classification string'

genomes = sorted(os.path.basename(d.rstrip('/')) for d in glob.glob(os.path.join(A, '*/')))
arch = {m for m in genomes if domain_of(m) == 'Archaea'}
bact = [m for m in genomes if m not in arch]
print(f"  domain assignment route: {route}")
print(f"  {len(bact)} bacterial, {len(arch)} archaeal genome(s)")

# ------------------------------------------------------- step 26 in-situ status
breadth = {}
cx = os.path.join(X, 'cross_abundance.tsv')
if os.path.exists(cx):
    byid = {mid(m): m for m in genomes}
    for r in csv.DictReader(open(cx), delimiter='\t'):
        k = byid.get(r.get('mag_id', ''))
        if k:
            try: breadth[k] = float(r.get('sed_breadth') or 0)
            except ValueError: pass
    print(f"  step 26 sediment breadth loaded for {len(breadth)}/{len(genomes)} genomes")
else:
    print("  step 26 cross_abundance.tsv not found — no in-situ weighting applied")
def insitu(m): return breadth.get(m, None)

rows = []   # one dict per genome

# =============================================================================
if MODE == 'antismash':
# =============================================================================
    reg = defaultdict(list)
    rp = os.path.join(B, 'bgc_regions.tsv')
    if os.path.exists(rp):
        for r in csv.DictReader(open(rp), delimiter='\t'):
            reg[r['genome']].append(r)
    def measured(m):
        d = os.path.join(B, 'as', m)
        return bool(glob.glob(os.path.join(d, '*.json'))) or \
               os.path.exists(os.path.join(d, 'index.html'))

    rule('1. antiSMASH REGIONS PER GENOME')
    nm = [m for m in genomes if not measured(m)]
    if nm:
        print(f"  !! {len(nm)} genome(s) have no antiSMASH output: "
              f"{', '.join(mid(x) for x in nm)}")
        print("     Their rows below are NOT MEASURED, not zero.")
        print()
    print(f"  {'genome':<30}{'dom':<5}{'sed.br':>8}{'regions':>9}  products")
    print('-' * W)
    for m in genomes:
        if not measured(m):
            rows.append(dict(genome=m, mag_id=mid(m), genus=genus(m),
                             domain=domain_of(m), sed_breadth=insitu(m),
                             tested=False, n_bgc=0, products=Counter()))
            print(f"  {mid(m):<30}{domain_of(m)[:4]:<5}{'':>8}   ---- not measured ----")
            continue
        rs = reg.get(m, [])
        prods = Counter()
        for r in rs:
            for p in (r['products'] or '').split(';'):
                if p: prods[p] += 1
        sb = insitu(m)
        rows.append(dict(genome=m, mag_id=mid(m), genus=genus(m),
                         domain=domain_of(m), sed_breadth=sb, tested=True,
                         n_bgc=len(rs), products=prods))
        print(f"  {mid(m):<30}{domain_of(m)[:4]:<5}"
              f"{('%.2f' % sb) if sb is not None else '   -':>8}{len(rs):>9}  "
              f"{', '.join(f'{k}x{v}' for k, v in prods.most_common(4))[:50]}")
    allprod = Counter()
    for r in rows:
        for k, v in r['products'].items(): allprod[k] += v
    print()
    print("  product classes across the catalogue:")
    for k, v in allprod.most_common():
        print(f"    {k:<34}{v:>4}")

# =============================================================================
else:   # Pfam floor
# =============================================================================
    sel = {}
    sp = os.path.join(B, 'hmm', 'selected.tsv')
    if os.path.exists(sp) and os.path.getsize(sp) > 0:
        for r in csv.DictReader(open(sp), delimiter='\t'):
            sel[r['model']] = (r['class'], r['confidence'], r['desc'])
    pats = [(r['class'], r['confidence'], re.compile(r['pattern'], re.I))
            for r in csv.DictReader(open(os.path.join(B, 'hmm', 'patterns.tsv')),
                                    delimiter='\t')]
    def classify(model, desc):
        if model in sel: return sel[model][0], sel[model][1]
        t = f"{model} {desc}"
        for cls, conf, rx in pats:
            if rx.search(t): return cls, conf
        return None, None

    hits_by_genome = {}
    for m in genomes:
        dt = os.path.join(B, 'hmm', f'{m}.domtbl')
        if not os.path.exists(dt): continue
        prot = defaultdict(lambda: {'len': 0, 'cls': set(), 'models': set()})
        for line in open(dt, errors='replace'):
            if line.startswith('#'): continue
            f_ = line.split()
            if len(f_) < 23: continue
            if HMODE == 'scan':
                model, pname, plen = f_[0], f_[3], f_[5]
                desc = ' '.join(f_[22:])
            else:
                pname, plen, model = f_[0], f_[2], f_[3]
                desc = sel.get(model, ('', '', ''))[2]
            cls, conf = classify(model, desc)
            if not cls: continue
            try: plen = int(plen)
            except ValueError: plen = 0
            e = prot[pname]
            e['len'] = max(e['len'], plen)
            e['cls'].add(cls); e['models'].add(model)
        hits_by_genome[m] = dict(prot)

    # --------------------------------------------------------------- CONTROLS
    rule('1. CONTROLS')
    print("  POSITIVE — fatty acid synthesis is not optional in a bacterium, so")
    print("  the carrier protein (ACP/PP-binding) and the AMP-binding superfamily")
    print("  must be found in nearly every BACTERIAL genome. Archaea are excluded:")
    print("  they build isoprenoid ether lipids and do not run bacterial FAS-II,")
    print("  so an absence there is expected biology, not a failed search.")
    print()
    tested_b = [m for m in bact if m in hits_by_genome]
    def has(m, cls): return any(cls in v['cls'] for v in hits_by_genome.get(m, {}).values())
    n_carrier = sum(1 for m in tested_b if has(m, 'CARRIER'))
    n_amp     = sum(1 for m in tested_b if has(m, 'NRPS_A'))
    n_ks      = sum(1 for m in tested_b if has(m, 'PKS_KS'))
    nb = len(tested_b)
    for lab, k in (('CARRIER (ACP / PP-binding)', n_carrier),
                   ('AMP-binding superfamily', n_amp),
                   ('ketoacyl synthase (FabB/F)', n_ks)):
        frac = (k / nb * 100) if nb else 0
        print(f"    {lab:<32}{k:>3}/{nb:<4}{frac:>6.1f}% of bacterial genomes")
    CONTROL_OK = nb > 0 and (n_carrier / nb >= 0.70) and (n_amp / nb >= 0.70)
    print()
    if CONTROL_OK:
        print("    [PASS] the search sees primary-metabolism machinery where it must be.")
        print("           A zero for secondary metabolism below is therefore a tested zero.")
    else:
        print("    [FAIL] the positive control is absent from bacterial genomes that")
        print("           certainly have it. The search, the library or the thresholds")
        print("           are wrong. NO VERDICT will be printed.")
    print()
    tested_a = [m for m in arch if m in hits_by_genome]
    na_c = sum(1 for m in tested_a if has(m, 'CARRIER'))
    print(f"  ARCHAEA (reported, not controlled): CARRIER in {na_c}/{len(tested_a)}")
    print("    Whatever this number is, it is not evidence about the search.")

    # ------------------------------------------------- the architecture filter
    rule('2. THE DISCRIMINATOR — ARCHITECTURE, NOT DOMAIN PRESENCE')
    print("  Bacterial FAS-II is dissociated: one domain per protein, ~300-450 aa.")
    print("  Type-I PKS and NRPS are megasynthases: several core domains in ONE ORF")
    print(f"  of {MEGA}-5000 aa. If the filter is doing anything, the KS-carrying")
    print("  proteins should fall into two separated length populations.")
    print()
    kslens = []
    for m, prot in hits_by_genome.items():
        for p, v in prot.items():
            if 'PKS_KS' in v['cls'] and v['len']: kslens.append(v['len'])
    kslens.sort()
    if kslens:
        import statistics as st
        bins = [(0, 500), (500, 1000), (1000, 2000), (2000, 10 ** 7)]
        print(f"  {len(kslens)} KS-carrying proteins, median {st.median(kslens):.0f} aa")
        for lo, hi in bins:
            k = sum(1 for x in kslens if lo <= x < hi)
            bar = '#' * min(60, k)
            hi_s = '+' if hi > 10 ** 6 else str(hi)
            print(f"    {lo:>5}-{hi_s:<6} {k:>4}  {bar}")
        big = sum(1 for x in kslens if x >= MEGA)
        print()
        if big == 0:
            print(f"  No KS protein reaches {MEGA} aa. Every one of them is dissociated")
            print("  FAS-II. That is a clean negative for type-I PKS, and it is the")
            print("  filter demonstrating that it separates the two populations rather")
            print("  than a threshold asserted and then trusted.")
        else:
            print(f"  {big} KS protein(s) at or above {MEGA} aa — megasynthase candidates.")
    else:
        print("  No KS-carrying proteins at all. See the control above before reading")
        print("  that as biology.")

    # --------------------------------------------------------------- the calls
    rule('3. CALLS PER GENOME')
    print("  NRPS      >=1 condensation domain. The C-domain is the discriminating")
    print("            one; AMP-binding alone is a fatty-acid-CoA ligase.")
    print(f"  PKS-I     one ORF with KS + (AT or carrier) and >= {MEGA} aa.")
    print("  HYBRID    C and KS in the same ORF.")
    print("  Other classes are single-domain calls and are marked LOW confidence.")
    print()
    LOWCONF = {'TERPENE', 'RIPP_YCAO', 'RIPP_SPASM', 'BETALACTONE'}
    print(f"  {'genome':<28}{'dom':<5}{'sed.br':>7}{'NRPS':>6}{'PKS':>5}{'HYB':>5}"
          f"{'RiPP':>6}{'bact':>6}{'sid':>5}{'terp':>6}{'phn':>5}")
    print('-' * W)
    for m in genomes:
        prot = hits_by_genome.get(m)
        if prot is None:
            print(f"  {mid(m):<28}{domain_of(m)[:4]:<5}{'':>7}   ---- not searched ----")
            rows.append(dict(genome=m, mag_id=mid(m), genus=genus(m),
                             domain=domain_of(m), sed_breadth=insitu(m),
                             tested=False, n_bgc=0, calls=Counter()))
            continue
        c = Counter()
        for p, v in prot.items():
            cl, L = v['cls'], v['len']
            if 'NRPS_C' in cl and 'PKS_KS' in cl and L >= MEGA: c['HYBRID'] += 1; continue
            if 'NRPS_C' in cl: c['NRPS'] += 1
            if 'PKS_KS' in cl and ('PKS_AT' in cl or 'CARRIER' in cl) and L >= MEGA:
                c['PKS'] += 1
            for k in ('RIPP_LANTI', 'RIPP_YCAO', 'RIPP_SPASM'):
                if k in cl: c['RiPP'] += 1; break
            if 'BACTERIOCIN' in cl: c['BACTERIOCIN'] += 1
            if 'SIDEROPHORE' in cl: c['SIDEROPHORE'] += 1
            if 'TERPENE' in cl: c['TERPENE'] += 1
            if 'PHOSPHONATE' in cl: c['PHOSPHONATE'] += 1
        sb = insitu(m)
        total = c['NRPS'] + c['PKS'] + c['HYBRID'] + c['RiPP'] + c['BACTERIOCIN'] \
                + c['SIDEROPHORE'] + c['PHOSPHONATE']
        rows.append(dict(genome=m, mag_id=mid(m), genus=genus(m),
                         domain=domain_of(m), sed_breadth=sb, tested=True,
                         n_bgc=total, calls=c))
        print(f"  {mid(m):<28}{domain_of(m)[:4]:<5}"
              f"{('%.2f' % sb) if sb is not None else '   -':>7}"
              f"{c['NRPS']:>6}{c['PKS']:>5}{c['HYBRID']:>5}{c['RiPP']:>6}"
              f"{c['BACTERIOCIN']:>6}{c['SIDEROPHORE']:>5}{c['TERPENE']:>6}"
              f"{c['PHOSPHONATE']:>5}")
    print()
    print("  'terp' is excluded from the BGC total: squalene/phytoene synthases are")
    print("  membrane-lipid (hopanoid, archaeal isoprenoid) enzymes as often as they")
    print("  are secondary metabolism, and this method cannot tell those apart.")

    # ------------------------------------------------- the methanogen negative
    rule('4. NEGATIVE CONTROL — THE METHANOGENS')
    meth = [r for r in rows if 'methano' in (r['genus'] or '').lower()]
    if meth:
        tot = sum(r['calls']['NRPS'] + r['calls']['HYBRID'] for r in meth if r['tested'])
        for r in meth:
            print(f"    {r['mag_id']:<30}NRPS {r['calls']['NRPS']:>3}   "
                  f"hybrid {r['calls']['HYBRID']:>3}   total BGC {r['n_bgc']:>3}")
        print()
        if tot == 0:
            print("    [PASS] no NRPS condensation domains in any methanogen, which is")
            print("           what a methanogen should look like. The filter is not")
            print("           calling secondary metabolism wherever it looks.")
        else:
            print(f"    [WARN] {tot} NRPS call(s) in methanogens. Not impossible —")
            print("           some archaea do carry NRPS-like enzymes — but check the")
            print("           models behind these before reporting them.")
    else:
        print("    No genus containing 'Methano' in the catalogue; control skipped.")
    globals()['CONTROL_OK'] = CONTROL_OK

# =============================================================================
# shared: in-situ weighting, table, verdict
# =============================================================================
rule('5. WHAT ANY OF THIS SAYS ABOUT THE LAGOON')
have_b = [r for r in rows if r.get('sed_breadth') is not None]
if not have_b:
    print("  No step-26 breadth available, so every genome is weighted equally.")
    print("  That overstates the lagoon relevance of the enrichment-only genomes.")
    ins, enr = rows, []
else:
    ins = [r for r in have_b if r['sed_breadth'] >= BINSITU]
    enr = [r for r in have_b if r['sed_breadth'] <  BINSITU]
    print(f"  {len(ins)} genome(s) at sediment breadth >= {BINSITU:.2f} — real residents.")
    print(f"  {len(enr)} genome(s) below it — enrichment-specific, and a BGC in one of")
    print("  them is a fact about a bottle, not about Etoliko.")
    print()
    for lab, grp in (('IN SITU', ins), ('ENRICHMENT-ONLY', enr)):
        gm = [r for r in grp if r.get('tested')]
        tb = sum(r['n_bgc'] for r in gm)
        car = sum(1 for r in gm if r['n_bgc'] > 0)
        miss = len(grp) - len(gm)
        note = f"   ({miss} not measured)" if miss else ''
        print(f"    {lab:<18}{tb:>4} cluster/call(s) across {car}/{len(gm)} "
              f"measured genome(s){note}")
    print()
    top = sorted(ins, key=lambda r: -r['n_bgc'])[:6]
    if top and top[0]['n_bgc'] > 0:
        print("  Richest in-situ genomes:")
        for r in top:
            if r['n_bgc'] == 0: continue
            print(f"    {r['mag_id']:<30}{r['genus'][:22]:<24}"
                  f"breadth {r['sed_breadth']:.2f}   {r['n_bgc']} call(s)")
    else:
        print("  No in-situ genome carries a single call.")

# ------------------------------------------------------------------- the table
# If the control failed the numbers below are all zero and all meaningless. The
# METHOD COLUMN says so, in the file itself — a reader who opens the TSV without
# reading VERDICT.txt must not be able to mistake a broken search for an absence.
METHOD = MODE if globals().get('CONTROL_OK', True) else MODE + '_CONTROL_FAILED'
out = os.path.join(B, 'BGC_SUMMARY.tsv')
fields = ['mag_id', 'genome', 'genus', 'domain', 'sed_breadth', 'in_situ',
          'method', 'measured', 'n_bgc']
extra = ['NRPS', 'PKS', 'HYBRID', 'RiPP', 'BACTERIOCIN', 'SIDEROPHORE',
         'TERPENE', 'PHOSPHONATE'] if MODE != 'antismash' else []
with open(out, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(fields + extra + (['products'] if MODE == 'antismash' else []))
    for r in rows:
        sb = r.get('sed_breadth')
        base = [r['mag_id'], r['genome'], r['genus'], r['domain'],
                ('%.4f' % sb) if sb is not None else '',
                ('yes' if (sb is not None and sb >= BINSITU) else
                 ('no' if sb is not None else '')),
                METHOD, ('yes' if r.get('tested') else 'NO'), r['n_bgc']]
        if MODE == 'antismash':
            w.writerow(base + [';'.join(f'{k}:{v}' for k, v in r['products'].items())])
        else:
            w.writerow(base + [r['calls'][k] for k in extra])
print()
print(f"  written: {out}")

# ----------------------------------------------------------------- the verdict
rule('6. VERDICT')
V = []
ok = globals().get('CONTROL_OK', True)
if MODE == 'antismash':
    ms = [r for r in rows if r.get('tested')]
    tot = sum(r['n_bgc'] for r in ms)
    car = sum(1 for r in ms if r['n_bgc'] > 0)
    V.append(f"antiSMASH: {tot} biosynthetic region(s) across {car}/{len(ms)} genomes.")
    if len(ms) < len(rows):
        V.append(f"{len(rows) - len(ms)} genome(s) produced no antiSMASH output and are")
        V.append("excluded — their zeros are missing measurements, not absences.")
    it = sum(r['n_bgc'] for r in ins if r.get('tested'))
    V.append(f"Of these, {it} sit in genomes confirmed in the sediment (breadth >= {BINSITU:.2f}).")
    V.append("Cluster boundaries and product classes are antiSMASH's rule-based calls;")
    V.append("MIBiG similarity is in the per-genome HTML and has not been summarised here.")
elif not ok:
    V.append("NO VERDICT. The positive control failed: bacterial genomes that must")
    V.append("carry ACP and AMP-binding did not. Every zero below the control line is")
    V.append("a broken search, not an absent pathway. Nothing here may be reported.")
else:
    tot = sum(r['n_bgc'] for r in rows)
    car = sum(1 for r in rows if r['n_bgc'] > 0)
    it  = sum(r['n_bgc'] for r in ins)
    V.append(f"Pfam core-domain FLOOR (antiSMASH was not available): {tot} call(s)")
    V.append(f"across {car}/{len(rows)} genomes; {it} of them in genomes confirmed in situ.")
    V.append("")
    V.append("This is a floor and must be reported as one. It counts biosynthetic")
    V.append("machinery, not clusters. It cannot give a boundary, cannot name a")
    V.append("product, and misses every RiPP whose only signature is a short")
    V.append("precursor peptide Bakta did not call. The controls establish that the")
    V.append("search works, so a zero is a tested zero — but a tested zero for CORE")
    V.append("DOMAINS, which is a weaker statement than a tested zero for clusters.")
    V.append("Running antiSMASH would strengthen, not overturn, this.")
for line in V: print('  ' + line)
open(os.path.join(B, 'VERDICT.txt'), 'w').write('\n'.join(V) + '\n')
print()
print(f"  written: {os.path.join(B, 'VERDICT.txt')}")
PY
RC=$?

# =============================================================================
# 4. VALIDATION
# =============================================================================
echo
echo "######## VALIDATION ########"
OK=0
check "interpretation completed" "$RC" || OK=1
check "BGC_SUMMARY.tsv written" $([ -s "$B/BGC_SUMMARY.tsv" ] && echo 0 || echo 1) || OK=1
check "VERDICT.txt written"     $([ -s "$B/VERDICT.txt" ]     && echo 0 || echo 1) || OK=1

NROW=$(($(wc -l < "$B/BGC_SUMMARY.tsv" 2>/dev/null || echo 1) - 1))
check "one row per genome ($NROW/$NGEN)" $([ "$NROW" -eq "$NGEN" ] && echo 0 || echo 1) || OK=1

# A row per genome is not the same as a measurement per genome. The last run
# passed this gate with 29 rows, every one of them unmeasured.
NMEAS=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="measured") c=i; next}
                    c && $c=="yes"' "$B/BGC_SUMMARY.tsv" 2>/dev/null | wc -l)
check "genomes actually measured ($NMEAS/$NGEN)" \
      $([ "$NMEAS" -gt 0 ] && echo 0 || echo 1) || OK=1
if [ "$NMEAS" -gt 0 ] && [ "$NMEAS" -lt "$NGEN" ]; then
  echo "  [warn] $((NGEN - NMEAS)) genome(s) unmeasured — their zeros are absences"
  echo "         of data, and the summary marks them measured=NO."
fi

if grep -q '^NO VERDICT' "$B/VERDICT.txt" 2>/dev/null; then
  echo "  [STOP] the control failed — nothing in $B may be reported."
  OK=1
fi

echo
if [ "$OK" -eq 0 ]; then
  echo "  33 complete (method: $MODE)."
  runlog "33_antismash complete — method=$MODE"
else
  echo "  33 finished with problems — read the FAILs above before using anything."
  runlog "**33_antismash finished with FAILs**"
fi
echo "  ended: $(date -Is)"
exit "$OK"
