# Kagome γ=1.272 source-assembly audit

## Result

The locked Kagome source assembly and the immutable exported MOF agree
exactly. The solver-free audit independently rebuilt and compared all 15,671
affine rows, all nine PSD blocks, all right-hand sides, and the objective.
There were zero coefficient mismatches.

Together with the strict exact/interval ray in
`KAGOME_EXACT_EXPORTED_MODEL_CERTIFICATE.md`, this establishes

```text
Δ_bulk ≤ 1.272
```

for the explicitly declared sign- and cyclic-spin-symmetry-restricted
infinite-volume KMS ground-state state-polynomial relaxation. It is not an
unrestricted-state result.

## Source and physical setup

The audit gates:

- SpectralGap base commit
  `a1171c906ff2cc2901e58c2426397a2f68c32bb7`;
- the committed patch SHA-256
  `332c0931ac810289aa3713af0948f259c01189270706af58b262d60d994d4abd`;
- SHA-256 identities of the four patched source files.

It reconstructs:

- `N=13`, `d=3`, `lso=5`, γ=`159/125`;
- triangles `(1,2,3)`, `(1,4,5)`, `(2,6,7)`, `(3,8,9)`,
  `(4,10,11)`, `(5,12,13)`;
- the first two triangles as the inner patch and no extra edges;
- `H = Σ_(i,j) Sᵢ·Sⱼ =
  ¼Σ_(i,j)(XᵢXⱼ + YᵢYⱼ + ZᵢZⱼ)`, comprising 54 exact
  coefficients `1/4`;
- removal of words with odd parity in any Pauli axis by the Kagome `isz`
  rule and cyclic `X→Y→Z` representatives from `reduce_perm`.

These restrictions define a symmetry-restricted infinite-volume KMS
ground-state class, not a finite periodic 13-site spectrum.

## Exact comparison

The audit independently assembles:

| component | dimensions/count |
|---|---:|
| state PSD bases | 271, 104 |
| gap PSD bases | 18, 17 |
| nine-site strengthening PSD blocks | 1, 9, 36, 84, 126 |
| stationarity variables | 20 |
| λ ordinal | 54,944 |
| total variables | 54,944 |
| affine rows compared | 15,671 |

For the five strengthening blocks, a separate exact computational-basis
Pauli action reconstructs the coefficient of every symmetric matrix
coordinate. This covers the source `posepsd9!` path instead of trusting only
the exported block inventory.

The intended-rational coefficient inventory is:

```text
-318/125, -2, -193/125, -159/125, -1, -68/125, -1/2, -1/4,
1/4, 1/2, 1, 159/125, 193/125, 2, 318/125, 4
```

All affine offsets are exactly zero. The objective is exactly maximize
`+λ`, with λ at ordinal 54,944. The canonical source-row SHA-256 is
`7de9f007983d79adcaa0997e7a8aea170b9f73bc050d32fa6731290f11176eb2`.

## Replay

From the repository root:

```bash
julia --project=julia-env --startup-file=no --history-file=no \
  tracks/polyopt/solutions/sdp-gap-seekers/scripts/\
audit_kagome_source_assembly.jl \
  .bohr-handoff/artifacts/certificate-audit-b1a1cad-20260728T102942Z/\
kagome-1.272/audit.mof.json.gz \
  .external/SpectralGap
```

The measured solver-free runtime was 16m45s on one local CPU, with about
3.1 GiB resident memory. No optimizer was constructed or invoked.

The transcript is:

```text
.bohr-handoff/kagome-source-assembly-audit.tsv
SHA-256 377dbc70028b6c97fcf6214346b9362d8ed8ac098365c04b81cc1fc3ca4f4c1b
```

Acceptance requires `rows_compared=15671`,
`coefficient_mismatches=0`, `objective=54944 => 1//1`,
`optimizer_invoked=false`, and `source_assembly_equal=true`.
