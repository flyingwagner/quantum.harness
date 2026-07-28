# Solver-free validation report

A Square J1-J2 conic feasibility model is now assembled and independently
parsed without optimization. The checks below cover geometry, exact Pauli
algebra, basis identity, exact `M/G/K` coefficients, canonical bytes, conic
rendering, exported-ray replay, and certificate post-processing. Solver status
by itself is never promoted to a physical gap claim.

## Julia unit suite

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/test/runtests.jl
```

Result on Julia 1.11.9: `671/671` checks passed without optimization.

```text
solver-free homogeneous conic-ray verifier   90
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
shared core canonical wire grammar           17
complex Hermitian to real PSD rendering        5
Square shared-core inventory declarations      9
solver-free Square conic render               20
exact core M/G/K pair algebra                23
Square J1-J2 core M/G/K source gate         125
small finite-patch ED construction oracle     3
solver-free status runner contract           58
```

The homogeneous-ray fixtures distinguish accepted improving rays from equality,
PSD, and objective-sign failures. They also cover affine recession directions,
scale invariance, the Kagome-size cancellation pathology, and fail-closed
validation of PSD-preserving diagonal-congruence variable maps.

## Square exact-core gate

`check_square_core_mgk.jl` exhausts the `L=1,d=2,g=1/2` structured-basis upper
triangles: 247,456 positive pairs and 28 gap pairs, producing the required
247,540 component records. Every lower-triangle coefficient is independently
recomputed from swapped inputs and equals the conjugate upper entry. The run
references 74,602 canonical scalar rows and invokes no solver.

The separate `H=Z`, basis `[X,Y]` fixture fixes the symmetrized commutator and
complex-packing signs. Float Hamiltonian coefficients and unapplied symmetry
metadata are rejected.

## Square conic preflight

The exact `L=1,d=2,g=1/2,gamma=1/10` source plan rendered in 6m38s to a
25 MiB MOF artifact with SHA-256
`8cc83c7ed497a940823b9b57433c3b6b50f2aaa20ae9330746fa0ba0ede60685`.
It contains 74,602 scalar variables, `L(1)=1`, three nonzero stationarity
equalities, complex PSD blocks 703 and 7 rendered as real PSD blocks 1406 and
14, and a feasibility objective. Independent MOI parsing recovered 74,602
variables, four equalities, and the two expected PSD dimensions.
Both emitter and reader report that optimization was not invoked.

A fresh exact source rebuild then replayed the full MOF: all 49 affine and
495,560 PSD coefficients, variable names/order, right-hand sides, cone
coordinates/dimensions, and feasibility objective matched. This is a complete
solver-free render audit, not just an inventory check.

## Canonical Square core inventory

The complete native `core_mgk` artifact contains 161,886,794 canonical bytes,
247,484 pair records, all 247,540 required component records, 247,824 nonzero
coefficients, and 74,602 referenced scalar rows. Its SHA-256 is
`1d466cea7256a71d03d50f84d58f2d088d69360bfeb5b18161108d10c16e8549`.
An independent two-minute validation passed byte-identical decode/re-encode,
the envelope and eight section hashes, all derived IDs and references,
Hermitian diagonal reality, exact nonzero/zero record rules, and complete
pair/component coverage. No optimizer was invoked.
Independent full builds at gamma `1/10` and `1/5` produced byte-for-byte
identical math artifacts and envelopes, proving the symbolic core identity is
separate from numerical threshold evaluation in this implementation.

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
- The follow-up uniform-PSD-block/row equilibration job `22988046` also
  returned `SLOW_PROGRESS` / unknown. Its back-transformed ray improved the
  original-model equality residual to `9.914902e-12`, but still failed the
  unchanged `1e-12` replay. This is numerical conditioning evidence, not an
  infeasibility result.
- Diagonal-congruence job `22988126` preserved PSD equivalence and reduced
  solve time/memory, but its back-transformed equality residual
  `7.049624e-11` was worse. It too remains rejected numerical evidence.

## Remaining boundary

- a conic-render envelope binding the exact source to the emitted MOF;
- a frozen legacy source mapping/event trace if Gate C compatibility with
  SpectralGap is required rather than a native-to-native tensor diff;
- a source-gated Square MOF/status runner and independent conic replay;
- any strictly audited Square infeasibility ray and resulting bulk-gap bound.

No current output supports the phrase “certified Square J1-J2 bulk-gap bound.”
