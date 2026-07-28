# Safe raw-status runner

This runner replaces the legacy scripts that compare a `certify_*` return value
to integer `1`/`0`. The checked SpectralGap patch returns a `NamedTuple`;
termination, primal, dual, and objective fields are preserved without assigning
physical meaning to `flag`.

Source gate: SpectralGap base
`a1171c906ff2cc2901e58c2426397a2f68c32bb7`, checked patch SHA-256
`332c0931ac810289aa3713af0948f259c01189270706af58b262d60d994d4abd`,
patched `src/SpectralGap.jl` SHA-256
`940cd72b9c4bea39b6daaefd2b9797c54df450cb0afbfb4d9322b2df1b3838bb`,
patched `src/sdp.jl` SHA-256
`4ea362723bd7601e67db3bc27f21a4a11506791ce2b9b82cb7e4a60b0f6bae10`.
The patch intentionally adds no package dependency: in particular it does not
import Clarabel, which is absent from SpectralGap's upstream `Project.toml`.

## Locked calculation

- TFIM: `H = -Σᵢ₌₁⁸ ZᵢZᵢ₊₁ + 0.5Σᵢ₌₁⁹ Xᵢ`, open boundary,
  `N=9`, `d=2`, `lso=6`, `γ ∈ {0.25,0.26}`.
- Kagome: six triangles
  `[[1,2,3],[1,4,5],[2,6,7],[3,8,9],[4,10,11],[5,12,13]]`,
  first two inner, no extra edges, `S=σ/2`, `N=13`, `d=3`, `lso=5`,
  `γ ∈ {1,1.2,1.26,1.28}`.

The Slurm array contains these six points only. There is no N=27 or d=4 path.
The requested resources are CPU-only: 32 CPU and 121.6 GB per cell, one cell
per array task.

## Local solver-free checks

```bash
julia --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/gap_status_runner.jl \
  --list-points

julia --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/gap_status_runner.jl \
  --model tfim --gamma 0.25 --dry-run

julia --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/test/runtests.jl

tracks/polyopt/solutions/sdp-gap-seekers/test/patch_replay.sh \
  /path/to/SpectralGap
```

`--dry-run` does not import SpectralGap and cannot construct a Mosek optimizer.
The fixture tests exercise both the old integer result and the new `NamedTuple`
result without loading SpectralGap or Mosek. The patch replay test checks the
exact upstream commit/tree, applies the patch in a disposable clone, verifies
all four source hashes, and confirms every package imported by
`src/SpectralGap.jl` is declared in its own `Project.toml`. Omit the optional
path to let the test clone the upstream repository.

## Independent exported-ray replay

An exported MOF model and its `audit.variables.tsv` ray can be checked without
constructing an optimizer or invoking Mosek:

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/verify_gap_ray.jl \
  /path/to/audit.mof.json.gz \
  /path/to/audit.variables.tsv
```

Exit status `0` means the scale-normalized affine-recession, cone, and
objective-direction checks all passed; `1` means at least one check failed.
The output names each failed check. `accepted_floating_point_ray` is an
independent floating-point replay, not an exact rational or interval proof.
The model/ray pair and verifier source revision must remain part of the audit
record.

## xH5 run sequence

Before real submission, confirm the Hamiltonians above and probe the current
partition/resource availability. Then, from the xH5 checkout:

```bash
RUN_ID=status-$(date -u +%Y%m%dT%H%M%SZ)
tracks/polyopt/solutions/sdp-gap-seekers/scripts/prepare_gap_status_run.sh \
  --run-id "$RUN_ID" --instantiate
```

Preparation is not scientific compute: it creates or verifies the bootstrap
`julia-env/Manifest.toml`, copies SpectralGap into a read-only per-run source
tree, and creates a separate per-run environment bound to that copy. It saves
the frozen `Project.toml`, `Manifest.toml`, `Pkg.status`, Julia binary
path/version/SHA-256, exact harness commit/tree, and source checksums. The six
cells execute only this frozen environment. If a Manifest is absent or any
source/runtime hash differs, preparation or the cell fails closed; an
environment without a Manifest is never called reproducible.
With `--instantiate`, the exact `.external/SpectralGap` checkout is registered
as a path dependency before `Pkg.instantiate()` generates the Manifest; override
it explicitly with `--spectralgap-path PATH` if the checkout is elsewhere.

From the local checkout, perform the required Slurm feasibility check via the
standard harness mechanism:

```bash
HARNESS_CLUSTER_PROFILE=scnet-xh5 \
scripts/harness_slurm.sh submit --test-only \
  --array 6 \
  --run-spec "results/gap-status/$RUN_ID/run_spec.json" \
  --command gap-status-array-contract-v1 \
  --script tracks/polyopt/solutions/sdp-gap-seekers/scripts/gap_status_array.sbatch \
  --partition xhacnormalb --time 01:00:00 --cpus 32
```

Only after the queue/resource choice and the exact test-only response have been
ratified, remove `--test-only` to submit. The command exports
`HARNESS_RUN_SPEC`; the dedicated wrapper derives the validated run ID from
that path.

After the array leaves `squeue`, save final accounting and checksums on xH5:

```bash
tracks/polyopt/solutions/sdp-gap-seekers/scripts/finalize_gap_status_run.sh \
  --run-id "$RUN_ID" --job-id "$JOB_ID"
```

Each `results/gap-status/$RUN_ID/cells/<cell>/result.jsonl` contains:

- model and complete point parameters;
- raw flag, termination, primal, dual, objective plus its availability, wall
  time, and exception with the complete stack trace;
- exact source, patch, environment, basis-input, and Hamiltonian fingerprints;
- explicit `unavailable` records for residual and witness, because the current
  API does not export them.

`sacct.provisional.txt` is captured inside each cell. Finalization requires six
terminal array-task records, six matching result/source gates, and matching
job/task IDs; it writes one explicit success/fail row per cell plus an aggregate
outcome. `sacct.final.txt` and `SHA256SUMS.final` are produced only after that
validation. Scheduler state alone is never treated as scientific evidence.
