#!/usr/bin/env bash
#SBATCH --partition=xhacnormalb
#SBATCH --job-name=gap-status
#SBATCH --array=1-6
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=3800M
#SBATCH --output=results/gap-status/slurm/%A_%a.out
#SBATCH --error=results/gap-status/slurm/%A_%a.err

# Exactly six independently auditable cells. There is deliberately no N=27 or
# d=4 branch here. Run prepare_gap_status_run.sh first; this job fails closed
# if the frozen Manifest or source gate is absent/mismatched.

set -euo pipefail
umask 077

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../../.." && pwd)"
cd "$REPO_ROOT"
JULIA_BIN="${JULIA_BIN:-$HOME/julia-1.11.5/bin/julia}"
RUN_ID="${GAP_STATUS_RUN_ID:-}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-}"

die() {
  echo "gap_status_array: $*" >&2
  exit 64
}

if [[ -z "$RUN_ID" && -n "${HARNESS_RUN_SPEC:-}" ]]; then
  RUN_SPEC_DIR="$(dirname "$HARNESS_RUN_SPEC")"
  if [[ "$(basename "$HARNESS_RUN_SPEC")" == "run_spec.json" &&
    "$(basename "$(dirname "$RUN_SPEC_DIR")")" == "gap-status" ]]; then
    RUN_ID="$(basename "$RUN_SPEC_DIR")"
  fi
fi

[[ "$RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,79}$ ]] ||
  die "set GAP_STATUS_RUN_ID or HARNESS_RUN_SPEC=results/gap-status/ID/run_spec.json"
[[ "$TASK_ID" =~ ^[1-6]$ ]] ||
  die "SLURM_ARRAY_TASK_ID must be an integer from 1 through 6"
[[ "${HARNESS_COMMAND:-}" == "gap-status-array-contract-v1" ]] ||
  die "HARNESS_COMMAND must be the audited gap-status-array-contract-v1 token"

RUN_DIR="$REPO_ROOT/results/gap-status/$RUN_ID"
RUN_SPEC="${HARNESS_RUN_SPEC:-$RUN_DIR/run_spec.json}"
[[ -s "$RUN_SPEC" ]] || die "missing run spec: $RUN_SPEC"

CELL_FIELDS="$(
  python3 - "$RUN_SPEC" "$RUN_ID" "$TASK_ID" <<'PY'
import json
import sys

path, run_id, task_text = sys.argv[1:]
expected = [
    ("01-tfim-gamma-0p25", "tfim", "0.25"),
    ("02-tfim-gamma-0p26", "tfim", "0.26"),
    ("03-kagome-gamma-1", "kagome", "1"),
    ("04-kagome-gamma-1p2", "kagome", "1.2"),
    ("05-kagome-gamma-1p26", "kagome", "1.26"),
    ("06-kagome-gamma-1p28", "kagome", "1.28"),
]
with open(path, encoding="utf-8") as handle:
    spec = json.load(handle)
if spec.get("schema_version") != "gap-status-run-spec-v1":
    raise SystemExit("run_spec schema mismatch")
if spec.get("run_id") != run_id:
    raise SystemExit("run_spec run_id mismatch")
cells = spec.get("cells")
actual = [
    (cell.get("cell_id"), cell.get("params", {}).get("model"), cell.get("params", {}).get("gamma"))
    for cell in cells
] if isinstance(cells, list) else None
if actual != expected:
    raise SystemExit("run_spec cells do not match the locked six-cell contract")
index = int(task_text) - 1
print("\t".join(expected[index]))
PY
)" || die "run spec validation failed"
IFS=$'\t' read -r CELL_ID MODEL GAMMA <<<"$CELL_FIELDS"
[[ -n "$CELL_ID" && -n "$MODEL" && -n "$GAMMA" ]] ||
  die "run spec selected an empty field"

if [[ "$MODEL" == "tfim" ]]; then
  CPUS_EXPECTED=16
elif [[ "$MODEL" == "kagome" ]]; then
  CPUS_EXPECTED=32
else
  die "run spec selected unsupported model: $MODEL"
fi

export HARNESS_RUN_SPEC="$RUN_SPEC"
FREEZE_DIR="$RUN_DIR/environment"
SOURCE_FREEZE_DIR="$RUN_DIR/source"
CELL_DIR="$RUN_DIR/cells/$CELL_ID"
RUNNER="$SCRIPT_DIR/gap_status_runner.jl"
RESULT="$CELL_DIR/result.jsonl"

[[ -s "$FREEZE_DIR/Project.toml" ]] ||
  die "missing frozen Project.toml; run prepare_gap_status_run.sh"
[[ -s "$FREEZE_DIR/Manifest.toml" ]] ||
  die "missing frozen Manifest.toml; unlocked dependencies are not accepted"
[[ -s "$RUN_DIR/source-preflight.json" ]] ||
  die "missing source preflight; run prepare_gap_status_run.sh"
[[ -s "$RUN_DIR/prepare.runmeta" ]] || die "missing prepare.runmeta"
[[ -s "$SOURCE_FREEZE_DIR/SHA256SUMS" ]] || die "missing frozen source checksums"

meta_value() {
  local key="$1"
  awk -F= -v key="$key" '$1 == key {sub(/^[^=]*=/, ""); print; found++}
    END {if (found != 1) exit 1}' "$RUN_DIR/prepare.runmeta"
}

EXPECTED_REPO_COMMIT="$(meta_value repo_commit)" ||
  die "invalid repo_commit in prepare.runmeta"
EXPECTED_REPO_TREE="$(meta_value repo_tree)" ||
  die "invalid repo_tree in prepare.runmeta"
EXPECTED_JULIA_BIN="$(meta_value julia_bin)" ||
  die "invalid julia_bin in prepare.runmeta"
EXPECTED_JULIA_VERSION="$(meta_value julia_version)" ||
  die "invalid julia_version in prepare.runmeta"
EXPECTED_JULIA_SHA256="$(meta_value julia_bin_sha256)" ||
  die "invalid julia_bin_sha256 in prepare.runmeta"
EXPECTED_RUN_SPEC_SHA256="$(meta_value run_spec_sha256)" ||
  die "invalid run_spec_sha256 in prepare.runmeta"

command -v "$JULIA_BIN" >/dev/null 2>&1 || die "Julia not found: $JULIA_BIN"
JULIA_BIN="$(readlink -f "$JULIA_BIN")"
export JULIA_LOAD_PATH="@:@stdlib"
JULIA_VERSION="$("$JULIA_BIN" --startup-file=no --history-file=no --version)"
JULIA_BIN_SHA256="$(sha256sum "$JULIA_BIN" | awk '{print $1}')"
[[ "$JULIA_BIN" == "$EXPECTED_JULIA_BIN" ]] ||
  die "Julia path differs from prepare: $JULIA_BIN"
[[ "$JULIA_VERSION" == "$EXPECTED_JULIA_VERSION" ]] ||
  die "Julia version differs from prepare: $JULIA_VERSION"
[[ "$JULIA_BIN_SHA256" == "$EXPECTED_JULIA_SHA256" ]] ||
  die "Julia binary hash differs from prepare"
[[ "$(git -C "$REPO_ROOT" rev-parse HEAD)" == "$EXPECTED_REPO_COMMIT" ]] ||
  die "harness commit differs from prepare"
[[ "$(git -C "$REPO_ROOT" rev-parse 'HEAD^{tree}')" == "$EXPECTED_REPO_TREE" ]] ||
  die "harness tree differs from prepare"
[[ -z "$(git -C "$REPO_ROOT" status --porcelain --untracked-files=no)" ]] ||
  die "tracked harness worktree became dirty after prepare"
[[ "$(sha256sum "$RUN_SPEC" | awk '{print $1}')" == "$EXPECTED_RUN_SPEC_SHA256" ]] ||
  die "run spec hash differs from prepare"
(
  cd "$FREEZE_DIR"
  sha256sum -c SHA256SUMS
) || die "frozen environment checksum failed"
(
  cd "$RUN_DIR/source/SpectralGap"
  sha256sum -c "$SOURCE_FREEZE_DIR/SHA256SUMS"
) || die "frozen SpectralGap source checksum failed"

export GAP_STATUS_EXPECTED_HARNESS_COMMIT="$EXPECTED_REPO_COMMIT"
export GAP_STATUS_EXPECTED_HARNESS_TREE="$EXPECTED_REPO_TREE"
export GAP_STATUS_EXPECTED_JULIA_BIN="$EXPECTED_JULIA_BIN"
export GAP_STATUS_EXPECTED_JULIA_VERSION="$EXPECTED_JULIA_VERSION"
export GAP_STATUS_EXPECTED_JULIA_SHA256="$EXPECTED_JULIA_SHA256"

mkdir -p "$CELL_DIR"
exec > >(tee "$CELL_DIR/stdout.log") \
  2> >(tee "$CELL_DIR/stderr.log" >&2)

{
  printf 'schema=gap-status-runmeta-v1\n'
  printf 'run_id=%s\n' "$RUN_ID"
  printf 'cell_id=%s\n' "$CELL_ID"
  printf 'model=%s\n' "$MODEL"
  printf 'gamma=%s\n' "$GAMMA"
  printf 'host=%s\n' "$(hostname)"
  printf 'started_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID:-unavailable}"
  printf 'slurm_array_job_id=%s\n' "${SLURM_ARRAY_JOB_ID:-unavailable}"
  printf 'slurm_array_task_id=%s\n' "$TASK_ID"
  printf 'slurm_cpus_per_task=%s\n' "${SLURM_CPUS_PER_TASK:-unavailable}"
  printf 'expected_min_cpus=%s\n' "$CPUS_EXPECTED"
  printf 'run_spec=%s\n' "$RUN_SPEC"
  printf 'run_spec_sha256=%s\n' "$(sha256sum "$RUN_SPEC" | awk '{print $1}')"
  printf 'harness_command=%s\n' "$HARNESS_COMMAND"
  printf 'repo_commit=%s\n' "$(git -C "$REPO_ROOT" rev-parse HEAD)"
  printf 'repo_tree=%s\n' "$(git -C "$REPO_ROOT" rev-parse 'HEAD^{tree}')"
  printf 'julia_bin=%s\n' "$JULIA_BIN"
  printf 'julia_version=%s\n' "$JULIA_VERSION"
  printf 'julia_bin_sha256=%s\n' "$JULIA_BIN_SHA256"
  printf 'project_sha256=%s\n' "$(sha256sum "$FREEZE_DIR/Project.toml" | awk '{print $1}')"
  printf 'manifest_sha256=%s\n' "$(sha256sum "$FREEZE_DIR/Manifest.toml" | awk '{print $1}')"
} >"$CELL_DIR/runmeta"

if [[ "${SLURM_CPUS_PER_TASK:-0}" =~ ^[0-9]+$ ]] &&
  ((SLURM_CPUS_PER_TASK < CPUS_EXPECTED)); then
  die "cell needs at least $CPUS_EXPECTED CPUs, got $SLURM_CPUS_PER_TASK"
fi

export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-$((CPUS_EXPECTED / 2))}"
export MOSEKBINDIR="${MOSEKBINDIR:-$HOME/mosek/mosek/11.2/tools/platform/linux64x86/bin}"

set +e
"$JULIA_BIN" --startup-file=no --history-file=no --project="$FREEZE_DIR" "$RUNNER" \
  --run-spec "$RUN_SPEC" --cell-index "$TASK_ID" \
  --repo-root "$REPO_ROOT" --environment-dir "$FREEZE_DIR" \
  --output "$RESULT"
RUNNER_RC=$?
set -e

printf 'runner_exit_code=%s\n' "$RUNNER_RC" >>"$CELL_DIR/runmeta"
printf 'finished_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"$CELL_DIR/runmeta"

# This in-job sacct snapshot is intentionally labelled provisional; the final
# accounting record must be collected after the array leaves squeue.
if command -v sacct >/dev/null 2>&1 && [[ -n "${SLURM_ARRAY_JOB_ID:-}" ]]; then
  sacct -j "${SLURM_ARRAY_JOB_ID}_${TASK_ID}" \
    --format=JobID,State,ExitCode,MaxRSS,Elapsed,AllocCPUS -P \
    >"$CELL_DIR/sacct.provisional.txt" 2>&1 || true
else
  printf 'unavailable: sacct command or Slurm job id missing\n' \
    >"$CELL_DIR/sacct.provisional.txt"
fi

(
  cd "$CELL_DIR"
  find . -maxdepth 1 -type f ! -name 'SHA256SUMS*' -printf '%P\0' |
    sort -z |
    xargs -0 -r sha256sum >SHA256SUMS.provisional
)

exit "$RUNNER_RC"
