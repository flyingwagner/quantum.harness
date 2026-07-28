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
`8c6106f`. The Kagome solver statuses remain floating-point numerical
evidence, but exact post-processing now proves a strict ray for the exported
MOF under intended rational coefficient reconstruction. The TFIM row combines
the same strict rational/interval standard with a completed exact
source-to-MOF coefficient audit.

| model | config | numerical transition | independent evidence | conclusion |
|---|---|---|---|---|
| 1D TFIM | `N=9, g=0.5, d=2, lso=6`, sign- and reflection-symmetric | `(0.25075,0.25125]` | γ=0.25125 ray projected onto 2,705 exact rational equalities; four PSD blocks proved by 256-bit directed interval LDLᵀ; all 2,705 source rows/objective match exactly | `Δ_bulk ≤ 0.25125` for the declared symmetry-restricted KMS state class |
| Kagome Heisenberg | `N=13, d=3, lso=5`, sign-symmetric | `(1.270,1.272]` | strict exact ray for intended-rational exported MOF: 15,671 exact rows, objective `119089//14841408527 > 0`, nine rigorous PSD blocks; source assembly not yet audited | exported-model certificate; physical `Δ_bulk ≤ 1.272` claim withheld pending source binding |

## Status (2026-07-28, certificate audit)

- **TFIM:** exact rational projection and rigorous PSD membership pass. The
  source audit reproduces the Hamiltonian support, sign/reflection state class,
  matrix orientation, 136 stationarity variables, every affine coefficient
  and right-hand side, four PSD blocks, and the `+λ` objective with zero
  mismatches.
- **Kagome:** the completed xH5 A/B confirms 4,887 exact duplicates
  materially change cost, scale, and equality residual while preserving the
  same stalled/unknown status; both rays fail replay. A separate uniform
  PSD-block/row scaling improved the original residual by about 3.3× but also
  stalled and failed replay at `9.9149e-12`. Diagonal congruence also stalled
  and regressed to `7.0496e-11`. Unit-improvement homogeneous feasibility
  removed arbitrary scale but failed at `9.9912e-11`; row equilibration changed
  it only to `9.1335e-11`. Exact post-processing subsequently produced a
  strict intended-rational exported-MOF ray. Do not state the physical
  `Δ_bulk ≤ 1.272` result until a source audit independently binds that MOF to
  the Kagome Hamiltonian, geometry, nine PSD blocks, and sign-symmetric state
  class.
- **Square J1-J2:** the exact `L=1,d=2,g=1/2,γ=1/10` source assembly,
  structured-basis `M/G/K` core, conic render, and canonical status envelope
  pass the solver-free audit. The envelope deliberately says `unsolved` and
  `optimizer_invoked=false`; no Square gap number exists.

## Open items

1. Build the solver-free Kagome source audit. It must independently reproduce
   all 15,671 affine rows, the objective, all nine PSD blocks, γ=`159/125`,
   patch geometry, Hamiltonian convention, and sign-symmetric state class.
2. Keep the Square status envelope distinct from a future immutable
   solver-result envelope. A serious Square solve requires separately
   ratified physical setup and witness replay.
