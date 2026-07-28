# Solver-free validation report

No Square J1-J2 conic model was assembled or solved. The checks below cover
geometry, exact Pauli algebra, basis identity, exact `M/G/K` pair coefficients,
exported-ray replay, and certificate post-processing. Solver status by itself
is never promoted to a physical gap claim.

## Julia unit suite

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/test/runtests.jl
```

Result on Julia 1.11.9: `589/589` checks passed without optimization.

```text
solver-free homogeneous conic-ray verifier   59
TFIM source-audit row comparison               6
square patch geometry                        24
status runner static safety gates            29
Pauli canonicalization                       10
bare Pauli basis counts                      72
full state-polynomial formal counts          13
storage estimates                             2
exact local spin identities                  22
generic solver-free problem adapter          43
structured basis manifests                  100
exact core M/G/K pair algebra                23
Square J1-J2 core M/G/K source gate         125
small finite-patch ED construction oracle     3
solver-free status runner contract           58
```

The homogeneous-ray fixtures distinguish accepted improving rays from equality,
PSD, and objective-sign failures. They also cover affine recession directions,
scale invariance, and the Kagome-size cancellation pathology.

## Square exact-core gate

`check_square_core_mgk.jl` exhausts the `L=1,d=2,g=1/2` structured-basis upper
triangles: 247,456 positive pairs and 28 gap pairs, producing the required
247,540 component records. Every lower-triangle coefficient is independently
recomputed from swapped inputs and equals the conjugate upper entry. The run
references 74,602 canonical scalar rows and invokes no solver.

The separate `H=Z`, basis `[X,Y]` fixture fixes the symmetrized commutator and
complex-packing signs. Float Hamiltonian coefficients and unapplied symmetry
metadata are rejected.

## Small ED oracle

For the finite 3×3 internal-bond Hamiltonian at `g=1/2`, two independent matrix
builders agree exactly. Hermiticity, trace, and total-`Sᶻ` commutator residuals
are zero; the ground residual is `3.34e-15`. Its first distinct finite-window
separation `0.6877583922161636` is only an algebra oracle and is not a bulk-gap
estimate.

## Certificate artifacts

- The supplied TFIM γ=0.25125 ray passes independent floating replay. Exact
  rational equality projection and directed 256-bit interval LDLᵀ prove a
  strict ray. The source audit independently reproduces all 2,705 exact affine
  rows, zero right-hand sides, block orientation, and `+λ` objective with zero
  mismatches. The resulting bound applies only to the declared sign- and
  reflection-symmetric KMS state class.
- The completed xH5 Kagome γ=1.272 original/dedup A/B returned
  `SLOW_PROGRESS` / unknown for both cells. Independent replay rejected both
  at `1e-12`, with normalized equality residuals `3.245719e-11` and
  `3.079840e-10`. Deduplication reduced time, memory, and scale but did not
  produce a certificate.

## Remaining boundary

- canonical shared-core byte records, IDs, envelope, and full tensor artifact;
- Square normalization, stationarity, affine right-hand sides, objective, and
  complex-to-real cone rendering;
- a source-gated Square MOF/status runner and independent conic replay;
- any strictly audited Square infeasibility ray and resulting bulk-gap bound.

No current output supports the phrase “certified Square J1-J2 bulk-gap bound.”
