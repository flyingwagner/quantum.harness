# Kagome γ=1.272 strict source-bound certificate

## Result

The immutable Kagome MOF at γ=1.272 now has a strict, optimizer-free conic
ray over the intended rational reconstruction of its floating coefficients.
The ray has:

- exactly zero residual on all 15,671 homogeneous affine rows;
- exact improving objective `119089//14841408527 > 0`;
- rigorous positive-semidefinite proofs for all nine PSD blocks.

The subsequent source audit independently reproduced every row, block,
objective coefficient, patch convention, and state restriction with zero
mismatches. The combined result establishes `Δ_bulk ≤ 1.272` for the
declared sign- and cyclic-spin-symmetry-restricted infinite-volume KMS
ground-state class. It is not an unrestricted-state claim.

## Immutable inputs and generated evidence

The input MOF is:

```text
.bohr-handoff/artifacts/certificate-audit-b1a1cad-20260728T102942Z/
kagome-1.272/audit.mof.json.gz
SHA-256 2e59c395f3af1cdaf762cd98af3507e3b4ee7b75e315af338d8bef5a8bbb3a53
```

The final intended-coefficient ray is:

```text
.bohr-handoff/kagome-1.272-source-rational-exact-ray.tsv
SHA-256 c85c4a1255ea9afe11b32b6ab01fafeb2e7da525a34678cb0427c90cd27f0fc1
size 1,534,857 bytes
```

The independent audit transcript is:

```text
.bohr-handoff/kagome-rationalization-20260728T1925Z/
independent-source-rational-ray-audit.tsv
SHA-256 9e61c904aed74ac939fe1885db538869d910232a0c072f34cdfcec203a10436a
```

Generated artifacts remain under `.bohr-handoff/` and are not committed.

## Construction

The supplied row-equilibrated floating ray was normalized by its maximum
absolute entry and rationalized with value tolerance `1e-20` and zero
threshold `8e-21`. Exact reverse peeling then gave a binary-coefficient
affine ray. Its third PSD block was rigorously indefinite, including a
negative exact order-three principal minor.

An exact affine-kernel correction was constructed from 150 free block-three
coordinates. Exact rational row reduction enforced six independent kernel
vectors with supports

```text
1,2,3;4,5,6;7,8,9;10,11,12;13,14,15;16,17,18
```

The resulting binary-coefficient ray has SHA-256
`52bbab9d2e1d21d2fef0c2f8fe8cf36d73849bdd5f1be8d785a611127698329a`.
Reconstructing MOF coefficients with tolerance `1e-14` initially exposed 46
nonzero exact residuals. Forty-six exact private-pivot corrections removed
all of them without changing the improving objective or the declared
block-three kernel.

## Rigorous PSD argument

Blocks 1, 2, and 4 pass 256-bit directed-interval LDLᵀ after exact zero rows
are removed. Blocks 5 through 9 are exactly zero. Their minimum certified
positive pivot lower bounds are:

| Block | Active dimension | Minimum pivot lower bound |
|---|---:|---:|
| 1 | 271 | `8.402943358234833e-6` |
| 2 | 91 | `8.625986647081080e-6` |
| 4 | 17 | `6.662163776682889e-6` |

For block 3, let `K` contain the six exact kernel vectors above and remove
coordinates `3,6,9,12,15,18`. The removed-coordinate submatrix of `K` is
nonsingular and `M K = 0` exactly. Therefore `M` is rationally congruent to
the retained 12×12 principal submatrix direct-summed with a six-dimensional
zero block. Directed-interval LDLᵀ proves the retained matrix positive
definite, with minimum pivot lower bound
`2.571142147119503e-6`. This proves the full 18×18 block positive
semidefinite.

No solver status, floating eigenvalue threshold, or relaxed verifier
tolerance enters this proof.

## Replay

From the repository root:

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/audit_exact_ray_psd.jl \
  .bohr-handoff/artifacts/certificate-audit-b1a1cad-20260728T102942Z/\
kagome-1.272/audit.mof.json.gz \
  .bohr-handoff/kagome-1.272-source-rational-exact-ray.tsv \
  3 '1,2,3;4,5,6;7,8,9;10,11,12;13,14,15;16,17,18' 1e-14
```

Acceptance requires:

- `nonzero_exact_residuals = 0`;
- `objective_improvement_positive = true`;
- every block reports `proved = true`;
- block 3 reports
  `exact_kernel_reduced_interval_ldlt_positive_semidefinite`;
- `optimizer_invoked = false`.

## Source binding

`KAGOME_SOURCE_AUDIT.md` records the completed solver-free source audit. It
reproduces all 15,671 rows, the exact +λ objective, γ=`159/125`, the
six-triangle/first-two-inner patch, all 54 Heisenberg coefficients, and all
nine PSD blocks—including the five `posepsd9!` strengthening blocks. Its
canonical source-row SHA-256 is
`7de9f007983d79adcaa0997e7a8aea170b9f73bc050d32fa6731290f11176eb2`;
the comparison reports zero mismatches and no optimizer invocation.
