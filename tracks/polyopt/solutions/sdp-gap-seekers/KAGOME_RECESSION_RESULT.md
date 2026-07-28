# Kagome γ=1.272 normalized-recession result

## Claim boundary

Unit-improvement homogeneous feasibility removed the unbounded ray scale but
did not meet the independent equality tolerance. The returned point was
rejected at normalized tolerance `1e-12`.

Kagome γ=1.272 remains numerical instability. Mosek's `FEASIBLE_POINT`
solution status is not an independently accepted ray, an infeasibility
certificate, or a certified bulk-gap bound.

## Reproducibility

- Slurm job: `22988185`, state `COMPLETED`, exit `0:0`.
- Source commit/tree:
  `b0c3e050599edfa175cc3d5a9bd64abaaa321990` /
  `2481899c7bdefc01b95bf87b0f868749d7d24ccd`.
- Immutable original MOF SHA-256:
  `2e59c395f3af1cdaf762cd98af3507e3b4ee7b75e315af338d8bef5a8bbb3a53`.
- Normalized-recession MOF SHA-256:
  `ce6ec074f7bee5e0a01b2b29b240758e945317c90709357469606304d28acceb`.
- Environment: Julia 1.11.5, MosekTools 0.15.10.
- Resources: one `xhacnormalb` node, 36 CPU, 128 GiB, one-hour limit.
- Accounting: 4m55s elapsed, batch MaxRSS 20,163,316 KiB.
- Evidence:
  `.bohr-handoff/xh5-results/kagome-recession-b0c3e05-20260728T172000Z/`.

Every source, environment, input, and generated-model hash gate passed. The
model retained 54,944 variables and nine PSD blocks and contained 15,672
equalities: 15,671 homogeneous source rows plus `cᵀx=1`.

## Solver and replay

Mosek terminated after 101.28s with `SLOW_PROGRESS`
(`MSK_RES_TRM_STALL`) and reported primal and dual `FEASIBLE_POINT`.
Independent replay against the immutable original MOF found:

| quantity | value | test |
|---|---:|---|
| normalized equality residual | `9.991205557433147e-11` | fail |
| normalized PSD violation | `1.1232546680640168e-21` | pass |
| objective direction | `1.0` | pass |
| normalized improving objective | `8.105338171473591e-6` | pass |
| ray scale | `123375.48154615676` | diagnostic |

A local replay reproduced equality and objective values exactly. Its
below-tolerance PSD violation was `5.132556935296638e-21`; the overall
equality-residual rejection was unchanged.

## Interpretation

The normalized problem accomplished its intended structural goal: the
candidate has finite scale and unit improvement rather than an arbitrary
`1e17–1e34` unbounded iterate. The remaining failure is the affine kernel.
High-precision diagnosis gives normalized equality residual
`9.991153581879089e-11` and maximum row-scaled residual
`2.4977883954697722e-11`, led by row 8,839 and the 17×17 block.

The smallest distinct follow-up is row equilibration applied to the normalized
recession feasibility model, without variable scaling or duplicate removal.
This preserves the same finite feasible set and isolates whether Mosek's row
accuracy can cross the independent gate.
