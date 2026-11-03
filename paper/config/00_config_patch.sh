# =============================================================================
# 00_config_patch.sh  --  APPEND THIS TO 00_config.sh
#
#   cat 00_config_patch.sh >> 00_config.sh
#
# WHY
#   A failed `conda install` into ${CONDA_ENVS_SHARED}/NB_Genome rolled
#   back and left conda's metadata inconsistent, so `conda run -n NB_Genome`
#   now fails with DirectoryNotACondaEnvironmentError. The environment's
#   binaries are intact -- only conda's bookkeeping is broken, and that
#   directory is not writable by this user so it cannot be repaired here.
#
#   These replacement functions stop using `conda run` altogether. They locate
#   the environment prefix directly and invoke the binary from its bin/ with
#   PATH and LD_LIBRARY_PATH set. This is how workspace_1 originally invoked
#   Flye (see flye.log: Cmd: ${CONDA_ENVS_SHARED}/NB_Genome/bin/flye).
#
#   Side benefits: faster (no conda run overhead), immune to any further conda
#   metadata damage, and works for environments owned by other users.
#
#   Bash uses the LAST definition of a function, so appending this file
#   overrides the earlier definitions without editing them.
# =============================================================================

# roots to search for conda environment prefixes, in order
export CONDA_ROOTS="${CONDA_ENVS_USER} ${CONDA_ENVS_SHARED}"

# resolve an env NAME to its prefix directory
env_prefix() {  # env_prefix <envname>
  local name="$1" root
  for root in $CONDA_ROOTS; do
    [ -d "$root/$name" ] && { echo "$root/$name"; return 0; }
  done
  return 1
}

# run a tool in whichever env the smoketest proved working, WITHOUT conda run
run_tool() {  # run_tool <tool> <args...>
  local tool="$1"; shift
  local var="ENV_${tool//-/_}"
  local env="${!var:-}"

  if [ -z "$env" ]; then
    echo "  !! no working env resolved for $tool" >&2; return 127
  fi

  # base PATH
  if [ "$env" = "-" ]; then
    command "$tool" "$@"; return $?
  fi

  local prefix
  prefix=$(env_prefix "$env") || { echo "  !! env prefix not found for '$env'" >&2; return 127; }

  if [ ! -x "$prefix/bin/$tool" ]; then
    echo "  !! $tool not executable at $prefix/bin/$tool" >&2; return 127
  fi

  PATH="$prefix/bin:$PATH" \
  LD_LIBRARY_PATH="$prefix/lib:${LD_LIBRARY_PATH:-}" \
  CONDA_PREFIX="$prefix" \
    "$prefix/bin/$tool" "$@"
}

# same approach for version capture
record_version() {  # record_version <label> <env|-> <command...>
  local label="$1"; shift; local env="$1"; shift
  local tool="$1"; shift
  local out prefix

  if [ "$env" = "-" ]; then
    out=$(command "$tool" "$@" 2>&1 | grep -viE 'warning|^$' | head -2 | tr '\n' ' ')
  else
    prefix=$(env_prefix "$env")
    if [ -n "$prefix" ] && [ -x "$prefix/bin/$tool" ]; then
      out=$(PATH="$prefix/bin:$PATH" LD_LIBRARY_PATH="$prefix/lib:${LD_LIBRARY_PATH:-}" \
            "$prefix/bin/$tool" "$@" 2>&1 | grep -viE 'warning|^$' | head -2 | tr '\n' ' ')
    else
      out="NOT FOUND in env '$env'"
    fi
  fi

  out=$(echo "$out" | sed 's/\t/ /g' | cut -c1-180)
  [ -f "$VERSIONS" ] || printf 'step\ttool\tversion\tenv\tdate\n' > "$VERSIONS"
  printf '%s\t%s\t%s\t%s\t%s\n' "${STEP:-NA}" "$label" "$out" "$env" "$(date -Is)" >> "$VERSIONS"
  echo "  [version] $label = $out"
}

export -f env_prefix run_tool record_version runlog check 2>/dev/null || true

echo "[config] patched: bypassing 'conda run', invoking binaries by env prefix"
