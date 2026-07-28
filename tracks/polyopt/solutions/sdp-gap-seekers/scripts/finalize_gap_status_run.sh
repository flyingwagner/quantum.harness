#!/usr/bin/env bash

# Run on the xH5 login node after the Slurm array has left squeue. Saves final
# sacct state and recomputes all artifact checksums. It does not submit, retry,
# cancel, or modify a solver result.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../../.." && pwd)"
RUN_ID=""
JOB_ID=""

die() {
  echo "finalize_gap_status_run: $*" >&2
  exit 64
}

while (($#)); do
  case "$1" in
    --run-id)
      (($# >= 2)) || die "--run-id requires a value"
      RUN_ID="$2"
      shift 2
      ;;
    --job-id)
      (($# >= 2)) || die "--job-id requires a value"
      JOB_ID="$2"
      shift 2
      ;;
    --help)
      echo "Usage: finalize_gap_status_run.sh --run-id ID --job-id INTEGER"
      exit 0
      ;;
    *)
      die "unknown argument: $1"
      ;;
  esac
done

[[ "$RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,79}$ ]] ||
  die "invalid --run-id"
[[ "$JOB_ID" =~ ^[0-9]+$ ]] || die "--job-id must be an integer"
command -v sacct >/dev/null 2>&1 || die "sacct is unavailable"
command -v squeue >/dev/null 2>&1 || die "squeue is unavailable"

RUN_DIR="$REPO_ROOT/results/gap-status/$RUN_ID"
[[ -d "$RUN_DIR/cells" ]] || die "missing cells directory: $RUN_DIR/cells"
[[ -s "$RUN_DIR/run_spec.json" ]] || die "missing run_spec.json"
[[ -s "$RUN_DIR/source-preflight.json" ]] || die "missing source-preflight.json"

if ! QUEUE_STATE="$(squeue -h -j "$JOB_ID" -o '%T' 2>/dev/null)"; then
  die "squeue failed for job $JOB_ID"
fi
[[ -z "$QUEUE_STATE" ]] ||
  die "job $JOB_ID is still present in squeue: $QUEUE_STATE"

SACCT_FAILURES=0
printf 'task_id\tjob_id\tstate\texit_code\toutcome\n' \
  >"$RUN_DIR/sacct-cell-outcomes.tsv"
sacct -j "$JOB_ID" \
  --format=JobID,JobName,Partition,State,ExitCode,MaxRSS,Elapsed,AllocCPUS -P \
  >"$RUN_DIR/sacct.final.txt"

for task_id in 1 2 3 4 5 6; do
  state="$(
    awk -F'|' -v id="${JOB_ID}_${task_id}" '$1 == id {print $4; exit}' \
      "$RUN_DIR/sacct.final.txt"
  )"
  [[ -n "$state" ]] || die "sacct has no array-task row for ${JOB_ID}_${task_id}"
  exit_code="$(
    awk -F'|' -v id="${JOB_ID}_${task_id}" '$1 == id {print $5; exit}' \
      "$RUN_DIR/sacct.final.txt"
  )"
  case "$state" in
    PENDING|RUNNING|CONFIGURING|COMPLETING|SUSPENDED|RESIZING|REQUEUED|REQUEUE_FED|SIGNALING|STAGE_OUT)
      die "array task ${JOB_ID}_${task_id} is non-terminal: $state"
      ;;
    COMPLETED)
      if [[ "$exit_code" == "0:0" ]]; then
        task_outcome=success
      else
        task_outcome=failed
        SACCT_FAILURES=$((SACCT_FAILURES + 1))
      fi
      ;;
    *)
      task_outcome=failed
      SACCT_FAILURES=$((SACCT_FAILURES + 1))
      ;;
  esac
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$task_id" "${JOB_ID}_${task_id}" "$state" "$exit_code" "$task_outcome" \
    >>"$RUN_DIR/sacct-cell-outcomes.tsv"
done

python3 - "$RUN_DIR" "$RUN_ID" "$JOB_ID" <<'PY' \
  >"$RUN_DIR/artifact-validation.txt"
import json
import pathlib
import sys

run_dir = pathlib.Path(sys.argv[1])
run_id = sys.argv[2]
job_id = sys.argv[3]
expected = [
    ("01-tfim-gamma-0p25", "tfim", "0.25"),
    ("02-tfim-gamma-0p26", "tfim", "0.26"),
    ("03-kagome-gamma-1", "kagome", "1"),
    ("04-kagome-gamma-1p2", "kagome", "1.2"),
    ("05-kagome-gamma-1p26", "kagome", "1.26"),
    ("06-kagome-gamma-1p28", "kagome", "1.28"),
]

spec = json.loads((run_dir / "run_spec.json").read_text(encoding="utf-8"))
actual_spec = [
    (cell.get("cell_id"), cell.get("params", {}).get("model"), cell.get("params", {}).get("gamma"))
    for cell in spec.get("cells", [])
]
if spec.get("schema_version") != "gap-status-run-spec-v1" or spec.get("run_id") != run_id:
    raise SystemExit("invalid run_spec header")
if actual_spec != expected:
    raise SystemExit("run_spec cells differ from locked contract")

preflight = json.loads((run_dir / "source-preflight.json").read_text(encoding="utf-8"))
if preflight.get("record_type") != "source_preflight":
    raise SystemExit("invalid source-preflight record type")
if preflight.get("source", {}).get("exact_source_match") is not True:
    raise SystemExit("source preflight did not pass")

cell_dirs = sorted(path.name for path in (run_dir / "cells").iterdir() if path.is_dir())
expected_dirs = sorted(cell_id for cell_id, _, _ in expected)
if cell_dirs != expected_dirs:
    raise SystemExit(f"cell directory set mismatch: {cell_dirs!r}")

exception_count = 0
nonzero_runner_count = 0
run_spec_sha256 = __import__("hashlib").sha256(
    (run_dir / "run_spec.json").read_bytes()
).hexdigest()
for cell_id, model, gamma in expected:
    cell_dir = run_dir / "cells" / cell_id
    metadata = {}
    for line in (cell_dir / "runmeta").read_text(encoding="utf-8").splitlines():
        key, separator, value = line.partition("=")
        if not separator or key in metadata:
            raise SystemExit(f"invalid runmeta in {cell_id}")
        metadata[key] = value
    if metadata.get("run_id") != run_id or metadata.get("cell_id") != cell_id:
        raise SystemExit(f"runmeta identity mismatch in {cell_id}")
    if metadata.get("slurm_array_job_id") != job_id:
        raise SystemExit(f"job id mismatch in {cell_id}")
    if metadata.get("slurm_array_task_id") != str(expected.index((cell_id, model, gamma)) + 1):
        raise SystemExit(f"array task id mismatch in {cell_id}")
    if metadata.get("run_spec_sha256") != run_spec_sha256:
        raise SystemExit(f"run spec hash mismatch in {cell_id}")
    if metadata.get("runner_exit_code") != "0":
        nonzero_runner_count += 1

    result_path = cell_dir / "result.jsonl"
    lines = result_path.read_text(encoding="utf-8").splitlines()
    if len(lines) != 1:
        raise SystemExit(f"result.jsonl must contain exactly one record in {cell_id}")
    result = json.loads(lines[0])
    required = {
        "model", "parameters", "flag", "termination", "primal", "dual",
        "objective", "objective_availability", "walltime", "exception", "source", "basis",
        "hamiltonian", "residual", "witness",
    }
    if result.get("schema_version") != "gap-status-result-v1":
        raise SystemExit(f"result schema mismatch in {cell_id}")
    if result.get("record_type") != "solver_result" or not required.issubset(result):
        raise SystemExit(f"incomplete result schema in {cell_id}")
    if result.get("model") != model or result.get("parameters", {}).get("gamma") != gamma:
        raise SystemExit(f"result parameters mismatch in {cell_id}")
    if result.get("source", {}).get("exact_source_match") is not True:
        raise SystemExit(f"source gate absent in {cell_id}")
    runtime = result.get("runtime", {})
    if runtime.get("slurm_array_job_id") != job_id:
        raise SystemExit(f"result job id mismatch in {cell_id}")
    if runtime.get("slurm_array_task_id") != str(expected.index((cell_id, model, gamma)) + 1):
        raise SystemExit(f"result task id mismatch in {cell_id}")
    if runtime.get("run_spec_sha256") != run_spec_sha256:
        raise SystemExit(f"result run spec hash mismatch in {cell_id}")
    for section in ("source", "basis", "hamiltonian"):
        if not result.get(section, {}).get("fingerprint"):
            raise SystemExit(f"missing {section} fingerprint in {cell_id}")
    if result.get("residual", {}).get("availability") != "unavailable":
        raise SystemExit(f"residual availability is not explicit in {cell_id}")
    if result.get("witness", {}).get("availability") != "unavailable":
        raise SystemExit(f"witness availability is not explicit in {cell_id}")
    if result.get("exception") is not None:
        exception_count += 1

print("validation=complete")
print(f"cells={len(expected)}")
print(f"exception_cells={exception_count}")
print(f"nonzero_runner_cells={nonzero_runner_count}")
PY

EXCEPTION_CELLS="$(
  awk -F= '$1 == "exception_cells" {print $2}' \
    "$RUN_DIR/artifact-validation.txt"
)"
NONZERO_RUNNER_CELLS="$(
  awk -F= '$1 == "nonzero_runner_cells" {print $2}' \
    "$RUN_DIR/artifact-validation.txt"
)"
[[ "$EXCEPTION_CELLS" =~ ^[0-9]+$ ]] ||
  die "artifact validation did not report exception count"
[[ "$NONZERO_RUNNER_CELLS" =~ ^[0-9]+$ ]] ||
  die "artifact validation did not report runner exit count"

FINAL_OUTCOME=success
if ((SACCT_FAILURES > 0 || EXCEPTION_CELLS > 0 || NONZERO_RUNNER_CELLS > 0)); then
  FINAL_OUTCOME=failed
fi
{
  printf 'schema=gap-status-finalize-v1\n'
  printf 'run_id=%s\n' "$RUN_ID"
  printf 'job_id=%s\n' "$JOB_ID"
  printf 'finalized_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'outcome=%s\n' "$FINAL_OUTCOME"
  printf 'sacct_failed_cells=%s\n' "$SACCT_FAILURES"
  printf 'exception_cells=%s\n' "$EXCEPTION_CELLS"
  printf 'nonzero_runner_cells=%s\n' "$NONZERO_RUNNER_CELLS"
} >"$RUN_DIR/finalize.runmeta"

find "$RUN_DIR/cells" -mindepth 1 -maxdepth 1 -type d -print0 |
  while IFS= read -r -d '' cell_dir; do
    (
      cd "$cell_dir"
      find . -maxdepth 1 -type f ! -name 'SHA256SUMS*' -printf '%P\0' |
        sort -z |
        xargs -0 -r sha256sum >SHA256SUMS
    )
  done

(
  cd "$RUN_DIR"
  find . -type f ! -name SHA256SUMS.final -printf '%P\0' |
    sort -z |
    xargs -0 -r sha256sum >SHA256SUMS.final
)

printf 'finalized %s outcome=%s\n' "$RUN_DIR" "$FINAL_OUTCOME"
[[ "$FINAL_OUTCOME" == "success" ]] || exit 1
