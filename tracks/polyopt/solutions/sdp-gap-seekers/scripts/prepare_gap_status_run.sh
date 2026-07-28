#!/usr/bin/env bash

# Prepare and freeze the exact Julia environment before submitting any status
# cells. This script never calls sbatch or a solver. `--instantiate` is explicit
# because it may download/resolve packages and generate julia-env/Manifest.toml.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../../.." && pwd)"
JULIA_BIN="${JULIA_BIN:-$HOME/julia-1.11.5/bin/julia}"
RUN_ID=""
INSTANTIATE=0
DRY_RUN=0
SPECTRALGAP_PATH=""

die() {
  echo "prepare_gap_status_run: $*" >&2
  exit 64
}

usage() {
  cat <<'EOF'
Usage:
  prepare_gap_status_run.sh --run-id ID [--instantiate]
                            [--spectralgap-path PATH] [--dry-run]

Without --instantiate, an existing julia-env/Manifest.toml is required.
With --instantiate, Pkg.develop(path=<SpectralGap>) plus Pkg.instantiate() may
generate/update that Manifest before the Project/Manifest/Pkg.status snapshot
is copied into results/. The default path is .external/SpectralGap.
No Slurm job and no Mosek optimization is started.
EOF
}

while (($#)); do
  case "$1" in
    --run-id)
      (($# >= 2)) || die "--run-id requires a value"
      RUN_ID="$2"
      shift 2
      ;;
    --instantiate)
      INSTANTIATE=1
      shift
      ;;
    --spectralgap-path)
      (($# >= 2)) || die "--spectralgap-path requires a value"
      SPECTRALGAP_PATH="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --help)
      usage
      exit 0
      ;;
    *)
      die "unknown argument: $1"
      ;;
  esac
done

[[ "$RUN_ID" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,79}$ ]] ||
  die "--run-id must match [A-Za-z0-9][A-Za-z0-9._-]{0,79}"

ENV_DIR="$REPO_ROOT/julia-env"
SPECTRALGAP_PATH="${SPECTRALGAP_PATH:-$REPO_ROOT/.external/SpectralGap}"
RUN_DIR="$REPO_ROOT/results/gap-status/$RUN_ID"
FREEZE_DIR="$RUN_DIR/environment"
SOURCE_FREEZE_DIR="$RUN_DIR/source"
FROZEN_SPECTRALGAP="$SOURCE_FREEZE_DIR/SpectralGap"
SLURM_LOG_DIR="$REPO_ROOT/results/gap-status/slurm"
RUNNER="$SCRIPT_DIR/gap_status_runner.jl"

[[ -f "$ENV_DIR/Project.toml" ]] || die "missing $ENV_DIR/Project.toml"

if ((DRY_RUN)); then
  printf 'PLAN run_id=%s\n' "$RUN_ID"
  printf 'PLAN repo=%s\n' "$REPO_ROOT"
  printf 'PLAN environment=%s\n' "$ENV_DIR"
  printf 'PLAN frozen_environment=%s\n' "$FREEZE_DIR"
  printf 'PLAN spectralgap_path=%s\n' "$SPECTRALGAP_PATH"
  printf 'PLAN frozen_spectralgap=%s\n' "$FROZEN_SPECTRALGAP"
  printf 'PLAN julia_bin=%s\n' "$JULIA_BIN"
  printf 'PLAN instantiate=%s\n' "$INSTANTIATE"
  printf 'PLAN freeze=%s\n' "$FREEZE_DIR"
  printf 'PLAN solver_called=false\n'
  exit 0
fi

command -v "$JULIA_BIN" >/dev/null 2>&1 || die "Julia not found: $JULIA_BIN"
JULIA_BIN="$(readlink -f "$JULIA_BIN")"
[[ -x "$JULIA_BIN" ]] || die "Julia is not executable: $JULIA_BIN"
export JULIA_LOAD_PATH="@:@stdlib"
JULIA_VERSION="$("$JULIA_BIN" --startup-file=no --history-file=no --version)"
JULIA_BIN_SHA256="$(sha256sum "$JULIA_BIN" | awk '{print $1}')"

TRACKED_DIRTY="$(git -C "$REPO_ROOT" status --porcelain --untracked-files=no)"
[[ -z "$TRACKED_DIRTY" ]] ||
  die "tracked harness worktree is dirty; commit the audited runner before preparing"
[[ -f "$SPECTRALGAP_PATH/Project.toml" ]] ||
  die "missing SpectralGap checkout: $SPECTRALGAP_PATH"
[[ -d "$SPECTRALGAP_PATH/.git" ]] ||
  die "SpectralGap must be a self-contained normal Git clone, not a linked worktree"

mkdir -p "$(dirname "$RUN_DIR")" "$SLURM_LOG_DIR"
mkdir "$RUN_DIR" ||
  die "run directory already exists; choose a fresh --run-id: $RUN_DIR"
mkdir "$FREEZE_DIR" "$SOURCE_FREEZE_DIR"

if ((INSTANTIATE)); then
  "$JULIA_BIN" --startup-file=no --history-file=no --project="$ENV_DIR" -e \
    'using Pkg; Pkg.develop(path=ARGS[1]); Pkg.instantiate(;verbose=true)' \
    "$SPECTRALGAP_PATH"
fi

[[ -f "$ENV_DIR/Manifest.toml" ]] ||
  die "Manifest.toml is missing; this environment is unlocked. Re-run with --instantiate, inspect the resolution, then freeze it."

cp "$ENV_DIR/Project.toml" "$FREEZE_DIR/Project.toml"
cp "$ENV_DIR/Manifest.toml" "$FREEZE_DIR/Manifest.toml"
cp -a "$SPECTRALGAP_PATH" "$FROZEN_SPECTRALGAP"
[[ -d "$FROZEN_SPECTRALGAP/.git" ]] ||
  die "frozen SpectralGap copy lost self-contained Git metadata"

# Rebind the copied Manifest to the immutable per-run source copy. This is the
# environment the cells execute; the mutable repo julia-env is never used by a
# scientific cell.
"$JULIA_BIN" --startup-file=no --history-file=no --project="$FREEZE_DIR" -e \
  'using Pkg; Pkg.develop(path=ARGS[1]); Pkg.instantiate(;verbose=true)' \
  "$FROZEN_SPECTRALGAP"
"$JULIA_BIN" --startup-file=no --history-file=no --project="$FREEZE_DIR" -e \
  'using Pkg; Pkg.status(;mode=Pkg.PKGMODE_MANIFEST)' \
  >"$FREEZE_DIR/Pkg.status.txt"

REPO_COMMIT="$(git -C "$REPO_ROOT" rev-parse HEAD)"
REPO_TREE="$(git -C "$REPO_ROOT" rev-parse 'HEAD^{tree}')"
export GAP_STATUS_EXPECTED_HARNESS_COMMIT="$REPO_COMMIT"
export GAP_STATUS_EXPECTED_HARNESS_TREE="$REPO_TREE"
export GAP_STATUS_EXPECTED_JULIA_BIN="$JULIA_BIN"
export GAP_STATUS_EXPECTED_JULIA_VERSION="$JULIA_VERSION"
export GAP_STATUS_EXPECTED_JULIA_SHA256="$JULIA_BIN_SHA256"
export MOSEKBINDIR="${MOSEKBINDIR:-$HOME/mosek/mosek/11.2/tools/platform/linux64x86/bin}"

{
  printf 'run_id=%s\n' "$RUN_ID"
  printf 'prepared_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'repo_commit=%s\n' "$REPO_COMMIT"
  printf 'repo_tree=%s\n' "$REPO_TREE"
  printf 'repo_dirty_tracked=0\n'
  printf 'julia_bin=%s\n' "$JULIA_BIN"
  printf 'julia_version=%s\n' "$JULIA_VERSION"
  printf 'julia_bin_sha256=%s\n' "$JULIA_BIN_SHA256"
  printf 'frozen_environment=%s\n' "$FREEZE_DIR"
  printf 'frozen_spectralgap=%s\n' "$FROZEN_SPECTRALGAP"
  printf 'lock_state=locked-by-manifest\n'
} >"$RUN_DIR/prepare.runmeta"

cat >"$RUN_DIR/run_spec.json" <<EOF
{
  "schema_version": "gap-status-run-spec-v1",
  "run_id": "$RUN_ID",
  "cells": [
    {"cell_id": "01-tfim-gamma-0p25", "params": {"model": "tfim", "gamma": "0.25"}},
    {"cell_id": "02-tfim-gamma-0p26", "params": {"model": "tfim", "gamma": "0.26"}},
    {"cell_id": "03-kagome-gamma-1", "params": {"model": "kagome", "gamma": "1"}},
    {"cell_id": "04-kagome-gamma-1p2", "params": {"model": "kagome", "gamma": "1.2"}},
    {"cell_id": "05-kagome-gamma-1p26", "params": {"model": "kagome", "gamma": "1.26"}},
    {"cell_id": "06-kagome-gamma-1p28", "params": {"model": "kagome", "gamma": "1.28"}}
  ]
}
EOF
printf 'run_spec_sha256=%s\n' \
  "$(sha256sum "$RUN_DIR/run_spec.json" | awk '{print $1}')" \
  >>"$RUN_DIR/prepare.runmeta"

(
  cd "$FREEZE_DIR"
  sha256sum Project.toml Manifest.toml Pkg.status.txt >SHA256SUMS
)
(
  cd "$FROZEN_SPECTRALGAP"
  sha256sum \
    Project.toml \
    src/SpectralGap.jl \
    src/basicfunction.jl \
    src/sdp.jl \
    src/strengthening.jl \
    >"$SOURCE_FREEZE_DIR/SHA256SUMS"
)

# Import and source/patch/Manifest verification only. The preflight path cannot
# construct an optimizer and therefore cannot consume a Mosek license.
"$JULIA_BIN" --startup-file=no --history-file=no --project="$FREEZE_DIR" "$RUNNER" \
  --run-spec "$RUN_DIR/run_spec.json" --cell-index 1 --preflight \
  --repo-root "$REPO_ROOT" --environment-dir "$FREEZE_DIR" \
  --output "$RUN_DIR/source-preflight.json"

chmod -R a-w "$FROZEN_SPECTRALGAP"
chmod -R a-w "$FREEZE_DIR"
printf 'prepared %s\n' "$RUN_DIR"
