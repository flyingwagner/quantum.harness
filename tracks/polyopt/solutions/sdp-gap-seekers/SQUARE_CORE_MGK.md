# Square J1-J2 exact core M/G/K boundary

This is a solver-independent source gate for the Square J1-J2 bulk-gap
hierarchy. It computes exact positive-moment `M`, commutator-energy `K`, and
covariance `G` coefficients for the materialized structured basis. It does not
construct a conic model, run a solver, or report a gap bound.

## Locked semantics

- Hamiltonian:
  `H = ¼Σ_J1(XX+YY+ZZ) + (g/4)Σ_J2(XX+YY+ZZ)`, with `J1=1` and exact rational
  `g`.
- Geometry: `Λ_L=[−L,L]²`, inner patch `I_L=Λ_(L−1)`. The outer patch is a
  local consistency window; no finite-volume boundary condition is imposed.
- State class: unrestricted infinite-volume KMS ground states, one flat block,
  and no symmetry quotient. A declared symmetry is rejected because the
  current manifest does not implement its action.
- Core convention:

  ```text
  M[j,k] = L(ζ(b_j† b_k))
  G[j,k] = L(ζ(a_j† a_k)) - L(ζ(a_j†)ζ(a_k))
  K[j,k] = 1/2 L(ζ(a_j†[H,a_k] - [H,a_j†]a_k))
  A_γ    = K - γG
  ```

`CoreMGK.jl` keeps `K`, `G_moment`, and `G_product` separate and derives the
symbolic γ coefficient of `A_γ`; numeric γ is not stored in the core tensor.
All input coefficients must be integers or rationals. Float Hamiltonians fail
closed.

The structured-manifest implementation is selectively carried from comparison
commit `7bacf012e3e775a95e4c042bc088d05d158cfc56`; the status-runner branch was
not merged wholesale.

## Solver-free source gate

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/check_square_core_mgk.jl
```

For the locked `L=1, d=2, g=1/2, γ=1/10` source fixture this checks every
upper-triangle pair and independently recomputes every swapped lower-triangle
coefficient. Current anchors are:

| quantity | exact/measured value |
|---|---:|
| Hamiltonian terms | 60 |
| positive basis | 703 rows |
| gap basis | 7 rows |
| positive upper pairs | 247,456 |
| gap upper pairs | 28 |
| required component records | 247,540 |
| nonzero row coefficients | 247,824 |
| referenced scalar rows | 74,602 |

Positive manifest SHA-256 is
`83befe24c09bccdc7d228fc60c606d301dd76c10688121e1e466d43a583d5c13`;
gap manifest SHA-256 is
`5be3d2db7be104d1bc431898496e8e34116787a7f14a30886fa6933924bea169`;
the source problem SHA-256 at γ=1/10 is
`f6f7cd7a0cc2e053e40ecd82f52a24438536869e3340b959cd7f68cab4467f4e`.

The independent hand fixture `H=Z`, basis `[X,Y]`, fixes the complex signs:
`M_XY=+iζ(Z)`, `K_XY=−2iζ(I)`, and
`G_XY=+iζ(Z)−ζ(X)ζ(Y)`. Upper-triangle `Tr(AQ)` packing therefore uses
`+2 Im(q)` for a `+i` coefficient, not a minus sign.

## Shortest path to a status/audit runner

The reviewed typed byte grammar, Pauli/basis/row content IDs, normalization,
stationarity, right-hand sides, feasibility objective, exact scalar-row
mapping, and complex-to-real PSD rendering are now implemented and tested.
The complete native canonical inventory and envelope now serialize all 247,540
component records and independently validate full coverage. The remaining
boundary is:

1. Add a conic-render envelope binding the core math hash, exact gamma,
   normalization/stationarity row IDs, real cone coordinates, MOF SHA, source
   commit, and environment.
2. Extend `gap_status_runner_lib.jl` with a Square point. The source
   gate must bind the Hamiltonian, `L`, `d`, basis manifests, core math hash,
   conic-render hash, state class, exact γ, environment, and output MOF. Raw
   solver status remains `unknown` unless an independently replayed witness or
   ray passes.
3. For a gating diff against legacy SpectralGap rather than another native
   emitter, freeze the explicit source mapping and complete source-event trace.

Acceptance for the eventual Square runner is: solver-free dry run reconstructs
every source hash; the emitted MOF independently reproduces its coefficient
tensor and exact setup; feasible status is never called a lower gap bound; and
an infeasible endpoint becomes `Δ_bulk ≤ γ` only after a strict witness audit.
Any future symmetry-restricted run needs an implemented, hashed action and a
separate state-class label; it cannot reuse the unrestricted result.
