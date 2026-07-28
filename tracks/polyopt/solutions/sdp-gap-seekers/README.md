# SDP Gap Seekers

## Team

| | |
|---|---|
| **Team name** | sdp-gap-seekers |
| **Members** | Xiansheng Cai (蔡贤盛), Sihan Hu (胡思寒) |

## Challenge

Certified bulk spectral-gap bounds for frustrated spin-1/2 models — compute
upper bounds on the locally non-degenerate bulk gap of infinite systems via the
state-polynomial SDP hierarchy of arXiv:2606.03836.

Addresses #88 — released by Xiangling Xu (许湘灵) and Jie Wang (王杰),
polyopt track.

## Approach

The first implementation target is the square-lattice J1-J2 Heisenberg model.
The finite patch is a local consistency window, not a periodic finite-volume
Hamiltonian. The hierarchy tests whether an infinite-volume KMS ground state
can have gap at least `gamma`; finite-level infeasibility excludes that
threshold, while finite-level feasibility does not prove gappedness.

Current sequence:

1. Specify the KMS/state-polynomial hierarchy and its result semantics.
2. Build deterministic square patches, Pauli-word reduction, and reproducible
   basis fingerprints.
3. Select a nested structured basis; a complete dense basis is already too
   large at small `(L,d)`.
4. Reproduce a published transverse-field Ising table entry with the same
   state class and symmetry restriction as the paper.
5. Assemble and validate Square J1-J2 at `g=0`, `0.50`, and `0.535`.
6. Add observable bounds and exact or interval certificate post-processing.

Shastry-Sutherland `g=0` remains the preferred positive-gap calibration after
the Square adapter is stable. Fallbacks #124 and #49 are not part of the current
implementation branch.

## Current artifacts

- [`square-j1j2-gap-sdp-spec.md`](square-j1j2-gap-sdp-spec.md): programmable
  mathematical specification and solver-status semantics.
- [`basis-counts.md`](basis-counts.md): exact solver-free formal basis counts
  and raw dense-memory estimates.
- [`local-identities.md`](local-identities.md): exact two-, three-, and
  four-site identities and their role in structured relaxations.
- [`spectralgap-refactor-plan.md`](spectralgap-refactor-plan.md): migration
  plan from model-specific code to a generic lattice/patch interface.
- [`validation-report.md`](validation-report.md): tests, finite-patch ED oracle,
  and the precise boundary of what has not yet been certified.
- [`STRICT_CERTIFICATE.md`](STRICT_CERTIFICATE.md): exact rational projection
  and interval-PSD post-processing for the supplied TFIM conic ray.
- [`TFIM_SOURCE_ASSEMBLY.md`](TFIM_SOURCE_ASSEMBLY.md): coefficient-exact
  source-to-MOF binding for the restricted TFIM relaxation.
- [`KAGOME_CONDITIONING.md`](KAGOME_CONDITIONING.md): solver-free γ=1.272
  conditioning diagnosis and exact equality-deduplication experiment.
- [`KAGOME_AB_RESULT.md`](KAGOME_AB_RESULT.md): completed xH5
  original-versus-deduplicated solve and unchanged-tolerance replay.
- [`structured-basis-manifest.md`](structured-basis-manifest.md): materialized,
  versioned Square basis rows and hashes.
- [`SQUARE_CORE_MGK.md`](SQUARE_CORE_MGK.md): exact Square `M/G/K` pair algebra,
  full solver-free pair-coverage gate, and the remaining conic-runner boundary.
- [`SHARED_CORE_WIRE.md`](SHARED_CORE_WIRE.md): canonical typed bytes and
  domain-separated content IDs shared by source and audit implementations.
- [`SQUARE_CONIC_RENDER.md`](SQUARE_CONIC_RENDER.md): exact normalization,
  stationarity, real PSD embedding, solver-free MOF command, and measured
  `L=1,d=2` inventory.
- [`SQUARE_CORE_INVENTORY.md`](SQUARE_CORE_INVENTORY.md): complete
  gamma-independent native `core_mgk` bytes, hash envelope, full
  pair/component coverage, and independent artifact validator.

The geometry/model layer uses exact arithmetic. The Square source gate can now
assemble and serialize a conic feasibility model through the declared JuMP/MOI
file API without invoking an optimizer. The external
Mosek/SpectralGap/QMBCertify environment reported in [`notes/`](notes/) is a
separate solver setup; its local patches must be committed and
regression-tested before this repository can rely on them.

## Result language

- `infeasible` is physically conclusive only with a valid solver status and
  auditable infeasibility evidence.
- timeout, numerical failure, and ambiguous statuses are `unknown`.
- floating-point results are numerical SDP bounds until certificate
  post-processing accounts for solver error.
- no current file reports a Square J1-J2 bulk-gap bound.

## Division of labor

- Solver environment, legacy SpectralGap fixes, and paper-baseline logs:
  Xiansheng.
- Square model specification, structured-basis prototype, exact algebra, and
  solver-independent tests: Sihan.
- SDP assembly, certificate validation, and reported scans: joint review.
