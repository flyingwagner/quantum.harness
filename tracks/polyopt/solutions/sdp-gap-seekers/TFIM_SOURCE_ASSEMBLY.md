# TFIM source-assembly equality audit

The locked TFIM export at `γ=0.25125 = 201/800` is coefficient-for-coefficient
equal to the patched SpectralGap source formulation. This audit closes the
source-assembly boundary left by `STRICT_CERTIFICATE.md`; it constructs no
optimizer and performs no solve.

## Reproduce

Use the patched SpectralGap tree at base commit
`a1171c906ff2cc2901e58c2426397a2f68c32bb7`:

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/audit_tfim_source_assembly.jl \
  /path/to/tfim-0.25125/audit.mof.json.gz \
  .external/SpectralGap
```

The command fails closed on a source/patch hash mismatch, any affine
coefficient or right-hand-side mismatch, a changed PSD block ordering, or a
changed objective. It reproduces the source rows with the same Pauli
reductions and compares them to the MOF after the explicit `1e-14` rational
coefficient reconstruction used by the strict post-processor.

## Locked formulation

- Hamiltonian:
  `H = -Σᵢ₌₁⁸ ZᵢZᵢ₊₁ + 1/2 Σᵢ₌₁⁹ Xᵢ`, open nine-site consistency window.
- Encoding: `3i-2 = Xᵢ`, `3i-1 = Yᵢ`, `3i = Zᵢ`.
- Hierarchy: `d=2`, `lso=6`, `γ=201/800`.
- State class: global-X-spin-flip invariant, so words with odd Y/Z parity
  vanish; open-chain reflection is also imposed by `reduce_mirror`.
- Semantics: this is a symmetry-restricted infinite-volume KMS ground-state
  state-polynomial relaxation, not a finite-chain excitation calculation.
- Matrix orientation: source upper triangles use `j≤k`; MOF coordinates use
  `PositiveSemidefiniteConeTriangle` `trimap`; the source applies the
  off-diagonal factor two exactly once.

## Exact comparison result

| field | result |
|---|---:|
| positive / gap basis dimensions | `211,50` / `11,14` |
| stationarity monomials | 136 |
| variables / homogeneous equalities | 23,949 / 2,705 |
| affine right-hand sides | all exact zero |
| objective | maximize `+λ`, ordinal 23,949 |
| coefficient rows compared | 2,705 |
| coefficient mismatches | 0 |
| canonical affine-row SHA-256 | `8ae45321ada0edbbc18f0761591cea6de6014e64c6e4a271f54181413b0aa88a` |

The exact coefficient inventory is

```text
-4, -2, -1, -201/400, -1/2, -201/800,
 201/800, 1/2, 201/400, 1, 2, 4.
```

Together with the exact affine projection and interval-PSD proof in
`STRICT_CERTIFICATE.md`, this proves infeasibility of this particular
symmetry- and reflection-restricted finite state-polynomial relaxation at
`γ=201/800`. Under the documented relaxation semantics it yields
`Δ_bulk ≤ 0.25125` for that restricted KMS state class. It is not an
unrestricted TFIM ground-state claim and not a statement about a pure-phase
quasiparticle gap.

