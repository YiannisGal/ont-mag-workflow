# SP07 shared helpers. Source with PART set to the part folder name.
export B=${DOCS_DIR}
export WS2=$B/Metagenome_analysis_PhD/workspace_2
export AMP=$B/Etoliko_138_2
export S=$WS2/99_thesis_audit/SP07_final
export ENVBIN=$S/env/.pixi/envs/default/bin
export CE=${CONDA_ENVS_USER}
export TMPDIR=/tmp/sp07_tmp; mkdir -p $TMPDIR   # short path: multiprocessing sockets fail (AF_UNIX path too long) under the SP07 path
export P=$S/$PART; cd $P
LOG=$P/COMMANDS.log
run() { echo "CMD [$(date -Is)] cwd=$PWD :: $*" >> $LOG; "$@"; local rc=$?; echo "EXIT [$(date -Is)] rc=$rc :: $1" >> $LOG; return $rc; }
md5in() { for f in "$@"; do printf '%s\t%s\n' "$(md5sum < "$f" | cut -c1-32)" "$f"; done >> $P/inputs_md5.tsv; }
ver() { printf '%s\t%s\t%s\n' "$1" "$2" "${3:-}" >> $P/versions.tsv; }
inenv() { local e=$1; shift; run env PATH=$CE/$e/bin:$PATH "$@"; }
