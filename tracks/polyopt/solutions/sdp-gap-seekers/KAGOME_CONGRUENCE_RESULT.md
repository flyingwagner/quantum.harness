# Kagome γ=1.272 diagonal-congruence result

Later exact affine/PSD post-processing and a zero-mismatch source audit
established `Δ_bulk ≤ 1.272` for the declared symmetry-restricted KMS class.
This page evaluates only the rejected floating congruence result.

## Claim boundary

Diagonal congruence balanced the large PSD blocks and reduced resource use,
but its back-transformed ray was rejected at the unchanged normalized
tolerance `1e-12`.

This floating result remains numerical instability. It is neither an
infeasibility certificate nor a certified upper bound on the bulk gap.

## Reproducibility

- Slurm job: `22988126`, state `COMPLETED`, exit `0:0`.
- Source commit/tree:
  `f05c95f4867213bb18e9edc9cf40b793d78c84a8` /
  `fcca3cc839a371cd03aeb10003d4a1fc8b09cf1d`.
- Immutable original MOF SHA-256:
  `2e59c395f3af1cdaf762cd98af3507e3b4ee7b75e315af338d8bef5a8bbb3a53`.
- Congruence MOF SHA-256:
  `4fe4358e695c7d850731508303a5bbb30996863f105cd651e3eb5e226caa6acc`.
- Scale-map SHA-256:
  `9217a927390029268a4e827c1a4397d17fbf2e41c2a1c6541876bc1b30fdec88`.
- Environment: Julia 1.11.5, MosekTools 0.15.10.
- Resources: one `xhacnormalb` node, 36 CPU, 128 GiB, one-hour limit.
- Accounting: 4m24s elapsed, batch MaxRSS 13,207,676 KiB.
- Evidence:
  `.bohr-handoff/xh5-results/kagome-congruence-f05c95f-20260728T170505Z/`.

All source, environment, input, generated-model, and scale-map hash gates
passed. The generated model retained 54,944 variables, 15,671 affine
equalities, and direct real PSD dimensions
`271,104,18,17,1,9,36,84,126`.

## Solver and replay

Mosek terminated after 45.00s with `SLOW_PROGRESS`
(`MSK_RES_TRM_STALL`), primal and dual `UNKNOWN_RESULT_STATUS`, and one
returned point.

The authoritative replay after exact back-transformation reported:

| quantity | value | test |
|---|---:|---|
| normalized equality residual | `7.049624122338534e-11` | fail |
| normalized PSD violation | `1.2730333745313615e-17` | pass |
| normalized improving objective | `1.6372512775795462e-6` | pass |
| original-coordinate scale | `6.66097141993352e17` | diagnostic |

An independent local replay reproduced the equality/objective values exactly
and reported normalized PSD violation `1.2621723971528551e-17`; the
below-tolerance eigenvalue difference did not affect the verdict.

## Comparison and decision

| formulation | equality residual | ray scale | solve wall | verdict |
|---|---:|---:|---:|---|
| original | `3.245719161e-11` | `4.319e17` | 90.93s | rejected |
| exact duplicate removal | `3.079839671e-10` | `3.619e16` | 31.11s | rejected |
| uniform PSD block + rows | `9.914901660e-12` | `5.808e34` | 66.56s | rejected |
| diagonal congruence + rows | `7.049624122e-11` | `6.661e17` | 45.00s | rejected |

The dominant residual returned to equality row 1 and the 271×271 block.
Balancing the 104×104 block did not stabilize the affine kernel. Repeating
reference-ray-derived coordinate maps is therefore low value.

The next structurally distinct experiment should normalize the homogeneous
recession problem itself: replace affine constraints by their recession
directions, impose improving objective `cᵀx=1`, and solve finite conic
feasibility. This is equivalent to existence of an improving homogeneous ray
but avoids asking the optimizer to represent an unbounded trajectory.
