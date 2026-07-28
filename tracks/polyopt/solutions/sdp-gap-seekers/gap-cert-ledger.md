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
`8c6106f`. These rows are floating-point evidence, not formal certificates.

| model | config | numerical transition | independent evidence | conclusion |
|---|---|---|---|---|
| 1D TFIM | `N=9, g=0.5, d=2, lso=6`, sign-symmetric | `(0.25075,0.25125]` | γ=0.25125 ray: equality `2.2785e-15`, PSD violation `1.1629e-22`, objective `7.0270e-6`, all normalized | strong replayable floating ray; no formal Δ bound yet |
| Kagome Heisenberg | `N=13, d=3, lso=5`, sign-symmetric | `(1.270,1.272]` | γ=1.272 ray rejected: normalized equality residual `6.6153e-11` at tolerance `1e-12`; variable scale `8.5896e16` | numerical instability; γ=1.272 is not infeasible or certified |

## Status (2026-07-28, certificate audit)

- **TFIM:** the exported ray is independently replayable at floating-point
  tolerance. Exact equalities plus rigorous PSD membership remain open.
- **Kagome:** the transition is numerical only. Do not move an upper bound
  through γ=1.272; the available ray fails the equality audit.
- **Square J1-J2:** no status/audit runner or gap number exists yet.

## Open items

1. Exact or interval post-process the accepted TFIM floating ray.
2. Run a source-locked Kagome conditioning A/B experiment; do not loosen the
   verifier tolerance.
3. Connect the structured Square basis to a source-gated coefficient assembly
   and three-way status/audit runner.
