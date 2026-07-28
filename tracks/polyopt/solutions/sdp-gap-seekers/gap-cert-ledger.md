# Gap-cert ledger — spectral-gap certificate audit

> The #88 challenge track: **certified upper bounds on the bulk spectral gap**
> Δ of (frustrated) spin-1/2 models, via the state-polynomial γ-feasibility SDP
> hierarchy of arXiv:2606.03836, using the turnkey `SpectralGap.jl` certifiers in
> `.external/SpectralGap` (`certify_Ising_gap`, `certify_Heisenberg_kagome_gap`).
>
> Direction (per the SPEC + arXiv:2606.03836): for a threshold γ, the SDP tests
> whether a KMS ground state can have a locally non-degenerate bulk gap ≥ γ.
> **Feasibility is monotone decreasing in γ**: small γ feasible, large γ
> infeasible. A strictly validated infeasible endpoint γ excludes a gap at
> least γ and gives Δ ≤ γ; the hierarchy converges downward as `(L,d)`
> increase. Orthogonality is encoded by the **covariance term**
> ω(a†a)−|ω(a)|², not an S=1 sector.
>
> Owner: xcai side. Kept separate from the energy-cert ledger
> (`feature/energy-cert-floor`) per the agreed energy/gap split.

## Methodology — γ-scan to locate the feasibility transition

The historical `certify_*(N, H, γ, d)` interface returned
`flag = (status==OPTIMAL ? 1 : 0)`:
- `flag=1` → γ feasible → Δ could be ≥ γ (not excluded).
- `flag=0` → some non-`OPTIMAL` outcome; no physical conclusion without the
  raw status and certificate audit.

A coarse γ-scan localizes the numerical transition. Only a validated
infeasibility witness at the upper endpoint can turn that endpoint into the
statement **Δ ≤ γ_upper**.

> ⚠️ The `flag=(status==OPTIMAL)` convention (SpectralGap.jl upstream) collapses
> all non-OPTIMAL statuses into flag=0. Per SPEC §8 this is unsafe for rigorous
> certification (timeouts/numerics also give flag=0). For run-time validation
> and bound-localization it is adequate; a residual/witness audit is needed
> before claiming a formally certified bound. (Todo: §8 audit.)

## Current audited ledger

The source solve is `b1a1cad`; the independent replay implementation is
`8c6106f`. The Kagome row remains floating-point evidence. The TFIM row now
combines strict rational/interval post-processing with an exact source-to-MOF
coefficient audit.

| model | config | numerical transition | independent evidence | conclusion |
|---|---|---|---|---|
| 1D TFIM | `N=9, g=0.5, d=2, lso=6`, sign- and reflection-symmetric | `(0.25075,0.25125]` | γ=0.25125 ray projected onto 2,705 exact rational equalities; four PSD blocks proved by 256-bit directed interval LDLᵀ; all 2,705 source rows/objective match exactly | `Δ_bulk ≤ 0.25125` for the declared symmetry-restricted KMS state class |
| Kagome Heisenberg | `N=13, d=3, lso=5`, sign-symmetric | `(1.270,1.272]` | six equivalent formulations rejected; best residual `9.9149e-12`, row-scaled unit-ray `9.1335e-11`, at tolerance `1e-12` | numerical instability; γ=1.272 is not infeasible or certified |

## Status (2026-07-28, certificate audit)

- **TFIM:** exact rational projection and rigorous PSD membership pass. The
  source audit reproduces the Hamiltonian support, sign/reflection state class,
  matrix orientation, 136 stationarity variables, every affine coefficient
  and right-hand side, four PSD blocks, and the `+λ` objective with zero
  mismatches.
- **Kagome:** the transition is numerical only. Do not move an upper bound
  through γ=1.272. The completed xH5 A/B confirms 4,887 exact duplicates
  materially change cost, scale, and equality residual while preserving the
  same stalled/unknown status; both rays fail replay. A separate uniform
  PSD-block/row scaling improved the original residual by about 3.3× but also
  stalled and failed replay at `9.9149e-12`. Diagonal congruence also stalled
  and regressed to `7.0496e-11`. Unit-improvement homogeneous feasibility
  removed arbitrary scale but failed at `9.9912e-11`; row equilibration changed
  it only to `9.1335e-11`.
- **Square J1-J2:** exact structured-basis `M/G/K` coefficient assembly now
  passes its full `L=1,d=2,g=1/2` solver-free pair/Hermiticity gate. No conic
  status/audit runner or gap number exists yet.

## Open items

1. Stop conditioning variants. Any further Kagome certificate attempt must
   exactly project the affine kernel and prove PSD preservation after
   correction; do not loosen the verifier tolerance.
2. Connect the structured Square basis to a source-gated coefficient assembly
   and three-way status/audit runner.
