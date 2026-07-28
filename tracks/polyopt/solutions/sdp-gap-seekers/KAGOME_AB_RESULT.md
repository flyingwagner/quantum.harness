# Kagome γ=1.272 exact-deduplication A/B

Later exact affine/PSD post-processing and a zero-mismatch source audit
established `Δ_bulk ≤ 1.272` for the declared symmetry-restricted KMS class.
This page remains the historical record of the rejected floating A/B rays.

The xH5 A/B confirms solver instability and does not certify infeasibility.
Both returned rays fail the unchanged independent normalized tolerance
`1e-12` because of affine equality residuals.

## Locked run

- Physical setup: spin-1/2 antiferromagnetic Kagome Heisenberg model on the
  N=13 six-triangle consistency patch, first two triangles inner, `d=3`,
  `lso=5`, sign-symmetric state class, `γ=1.272`.
- Source commit/tree:
  `729b3959e78b822191ce807b8cb50a1245e43179` /
  `aa09179833a43847b1b87e1810c3a3c78090e28b`.
- Original MOF SHA-256:
  `2e59c395f3af1cdaf762cd98af3507e3b4ee7b75e315af338d8bef5a8bbb3a53`.
- Deduplicated MOF SHA-256:
  `e67e8a4dfe66aaccd0e69917546aefb99e8ee68bf5db3cc038218a21280cedfc`.
- xH5 job: array `22986235`, `xhacnormalb`, account `giggleliu`, QoS
  `user_sihhu`, 36 CPUs and 128 GiB per cell, one-hour limit.
- Runtime: Julia 1.11.5, MosekTools 0.15.10, MOSEK license from the private
  xH5 account. Both tasks finished `COMPLETED/0:0`.

The initial 32-CPU submission was rejected before a job ID because 128 GiB
exceeded the scheduler's 3931-MiB-per-CPU rule. Array `22986209` then exercised
only the clean-tree source gate and failed in 1–2 seconds; it did not parse or
solve a model. Both failed attempts are retained in the handoff work log.

## Result

| quantity | original 15,671 rows | exact dedup 10,784 rows |
|---|---:|---:|
| MOSEK termination / raw status | `SLOW_PROGRESS` / `MSK_RES_TRM_STALL` | same |
| primal / dual result status | `UNKNOWN_RESULT_STATUS` / same | same |
| solve wall time | 90.933 s | 31.109 s |
| Slurm batch MaxRSS | 20,091,096 KiB | 7,044,028 KiB |
| ray max norm | `4.319157242331815e17` | `3.6194330301836056e16` |
| normalized equality residual | `3.2457191614730715e-11` | `3.0798396708916705e-10` |
| normalized PSD violation (xH5 replay) | `3.791858220170951e-21` | `4.225909493793811e-21` |
| normalized improving objective | `8.089662951328586e-6` | `8.089309335753688e-6` |
| verifier at `1e-12` | rejected: equality | rejected: equality |

Exact duplicate removal therefore reduces solve time by about 2.9×, peak
resident memory by about 2.85×, and ray scale by about 11.9×. It nevertheless
makes the scale-normalized equality defect about 9.5× worse. The nearly
unchanged normalized improving objective and tiny relative PSD violations,
combined with materially different equality residuals for exactly the same
feasible set, are direct evidence of formulation-sensitive numerical
instability.

Neither ray is eligible for the exact/interval post-processor: that step is
only entered after the floating verifier accepts. No Kagome infeasibility or
bulk-gap upper bound follows.

A second replay on the Bohrium host reproduced every decisive quantity and
both equality rejections. Its large-matrix eigensolver reported
`3.0852949622765852e-21` rather than `3.791858220170951e-21` for the original
relative PSD violation; both are many orders below tolerance and both pass the
cone check. The deduplicated PSD value was identical.

The 256-bit row diagnostic also reproduces the equality failure independently
of Float64 summation. Row 1 is worst in both formulations and is dominated by
variable 42,488 in the 17×17 second gap block. Its high-precision normalized
residual changes from `3.245619078841464e-11` to
`3.079816863624107e-10`, whereas Float64 summation differences are only
`432.27` and `82.55` in unnormalized units.

## Smallest next A/B

Keep the same locked physical setup, exact deduplicated row set, solver
tolerance, and verifier tolerance. Compare identity scaling against one
predeclared positive diagonal congruence/block equilibration map:

1. derive block scales from source coefficient norms, not from either failed
   ray;
2. record the invertible variable map and inverse coefficient map exactly;
3. prove the transformed MOF is algebraically equivalent by round-tripping
   every affine row, cone block, and objective;
4. solve both cells and replay mapped-back rays at `1e-12`.

Do not combine this with another row transformation or tolerance change; the
conditioning effect must remain identifiable.
