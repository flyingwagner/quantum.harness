# Strict TFIM conic-certificate post-processing

This post-processor turns the supplied TFIM floating ray at `γ=0.25125` into
an exact rational recession direction and proves its PSD conditions with
directed interval arithmetic. It performs no optimization and does not load
Mosek.

## Reproduce

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/postprocess_gap_ray.jl \
  /path/to/tfim-0.25125/audit.mof.json.gz \
  /path/to/tfim-0.25125/audit.variables.tsv \
  /path/to/tfim-0.25125-exact-ray.tsv
```

The command fails closed unless the model has 23,949 variables, 2,705 affine
equalities, PSD block dimensions `211,50,11,14`, and exactly the reconstructed
coefficient inventory

```text
-4, -2, -1, -201/400, -1/2, -201/800,
201/800, 1/2, 201/400, 1, 2, 4.
```

Here `γ=201/800` and `2γ=201/400`. The `1e-14` coefficient reconstruction
tolerance maps only the MOF's Float64 spellings of this displayed inventory;
the inventory assertion prevents silent reuse on a different model.

## Construction and checks

1. Scale the floating ray so `maxᵢ |xᵢ|=1`.
2. Reconstruct its coordinates as rationals within `1e-12`; entries at most
   `1e-20` in magnitude are set to exact zero.
3. Evaluate every rational affine equality. For each nonzero residual, select
   a pivot coordinate occurring in exactly that one equality and correct it
   exactly. Thus later corrections cannot disturb earlier rows.
4. Require all 2,705 corrected residuals to be exactly zero and the rational
   maximizing-objective direction to be strictly positive.
5. Reconstruct every symmetric PSD block over the rationals. Remove only rows
   that are identically zero over the rationals, then use a 256-bit,
   outward-rounded interval LDLᵀ factorization. A block passes only when every
   interval pivot has a strictly positive lower endpoint.

For the supplied artifact, 922 private-pivot corrections close all 922 initial
exact residuals. All four PSD blocks pass. The `11×11` block has exact zero
rows 8–11 and a rigorously positive-definite `7×7` principal block. The other
three blocks are rigorously positive definite. The corrected rational ray is
therefore a strict certificate for the explicitly reconstructed rational conic
model.

## Source binding and scope

`TFIM_SOURCE_ASSEMBLY.md` now closes the source-assembly gate. Its solver-free
audit recomputes all 2,705 affine rows from the source Hamiltonian, Pauli
reduction, sign/reflection restrictions, covariance orientation, stationarity
monomials, and `+λ` maximizing objective. Every reconstructed rational
coefficient and zero right-hand side matches the exported MOF; the canonical
row SHA-256 is
`8ae45321ada0edbbc18f0761591cea6de6014e64c6e4a271f54181413b0aa88a`.

The combined result is a strict infeasibility certificate for the locked
`N=9,d=2,lso=6,γ=201/800` sign- and reflection-symmetric state-polynomial
relaxation. It implies `Δ_bulk ≤ 0.25125` in that restricted infinite-volume
KMS ground-state class. It is not an unrestricted-state TFIM statement and is
not the finite/open-chain or pure-phase quasiparticle gap.
