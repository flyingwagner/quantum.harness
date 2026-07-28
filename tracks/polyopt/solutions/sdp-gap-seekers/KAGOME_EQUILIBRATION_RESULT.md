# Kagome γ=1.272 uniform-block equilibration result

## Claim boundary

The PSD-safe uniform-block/row conditioning experiment improved the equality
residual relative to the original solve, but its back-transformed ray was
rejected at the unchanged normalized tolerance `1e-12`.

Kagome γ=1.272 remains numerical instability. This result is neither an
infeasibility certificate nor a certified upper bound on the bulk gap.

## Reproducibility

- Slurm job: `22988046`, state `COMPLETED`, exit `0:0`.
- Source commit/tree:
  `57befde49ffa374bbadc15f032237305bb049a81` /
  `ae3e893c8668609e2552571b7091bfa7387f0041`.
- Immutable original MOF SHA-256:
  `2e59c395f3af1cdaf762cd98af3507e3b4ee7b75e315af338d8bef5a8bbb3a53`.
- Equilibrated MOF SHA-256:
  `830de74414ea3f1f16037376a9dfb3ec38fa01e31e5e80420f34c981a6411a83`.
- Scale-map SHA-256:
  `2ebfb643cf2b0ede6435e42686e8adb6de487963432f23043e555bf4cde68e4a`.
- Environment: Julia 1.11.5, MosekTools 0.15.10.
- Resources: one `xhacnormalb` node, 36 CPU, 128 GiB, one-hour limit.
- Accounting: 4m53s elapsed, batch MaxRSS 15,957,768 KiB.
- Evidence:
  `.bohr-handoff/xh5-results/kagome-equilibrated-57befde-20260728T162500Z/`.

The source, tree, scripts, environment, immutable input, generated model, and
scale map were all hash-gated before the solve. The generated model preflight
reported 54,944 variables, 15,671 affine equalities, and direct real PSD
dimensions `271,104,18,17,1,9,36,84,126`.

## Solver and replay

Mosek terminated after 66.56s with:

- termination `SLOW_PROGRESS` (`MSK_RES_TRM_STALL`);
- primal and dual `UNKNOWN_RESULT_STATUS`;
- one returned point.

The point was back-transformed and replayed against the immutable original
MOF. The authoritative xH5 replay reported:

| quantity | value | test |
|---|---:|---|
| normalized equality residual | `9.914901660490167e-12` | fail |
| normalized PSD violation | `1.2718060369584767e-17` | pass |
| normalized improving objective | `6.452369852106409e-8` | pass |
| original-coordinate scale | `5.808433831135018e34` | diagnostic |

A second local replay reproduced the equality and objective values exactly.
Its normalized PSD violation was `1.646984933093529e-17`, reflecting
platform/LAPACK sensitivity in an eigenvalue many orders below tolerance; its
PSD pass and overall equality-residual rejection were unchanged.

## Comparison

| formulation | equality residual | ray scale | solve wall | verdict |
|---|---:|---:|---:|---|
| original | `3.245719161e-11` | `4.319e17` | 90.93s | rejected |
| exact duplicate removal | `3.079839671e-10` | `3.619e16` | 31.11s | rejected |
| uniform PSD block + row scaling | `9.914901660e-12` | `5.808e34` | 66.56s | rejected |

The uniform map improved the original equality residual by about 3.3×, but
missed acceptance by about 9.9× and greatly increased the back-transformed
ray scale. Its largest residual rows remain cancellation-dominated by the
271×271 and 17×17 PSD blocks. Loosening the verifier tolerance is not a valid
response.

## Next experiment

The smallest distinct experiment is diagonal congruence scaling inside each
direct PSD block, `X=DZD`, followed by the same row/objective normalization.
It is exactly reversible, preserves PSD membership in both directions, and
can balance the 104×104 block's roughly 20-order diagonal range. Acceptance
remains a back-transformed replay against the immutable original MOF at
`1e-12`.
