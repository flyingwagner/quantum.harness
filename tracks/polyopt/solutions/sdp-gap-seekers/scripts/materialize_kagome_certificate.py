#!/usr/bin/env python3

"""Build or verify the fail-closed Kagome γ=1.272 certificate envelope."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path


PROOF_COMMIT = "0b58fc75dddd4830ae7b6b18a91bf61adff09d12"
PROOF_TREE = "3196131c33f04e1bcd8e8f1fb3f7f27785da90b4"
SPECTRALGAP_COMMIT = "a1171c906ff2cc2901e58c2426397a2f68c32bb7"

EXPECTED_ARTIFACTS = {
    "model": "2e59c395f3af1cdaf762cd98af3507e3b4ee7b75e315af338d8bef5a8bbb3a53",
    "ray": "c85c4a1255ea9afe11b32b6ab01fafeb2e7da525a34678cb0427c90cd27f0fc1",
    "exact_audit": "9e61c904aed74ac939fe1885db538869d910232a0c072f34cdfcec203a10436a",
    "source_audit": "377dbc70028b6c97fcf6214346b9362d8ed8ac098365c04b81cc1fc3ca4f4c1b",
}

PROOF_FILES = {
    "tracks/polyopt/solutions/sdp-gap-seekers/src/GapRayPostprocess.jl":
        "3ca59e6312db993eee768e6566b28f41865456eae8b84f92f1edda9c13d509fa",
    "tracks/polyopt/solutions/sdp-gap-seekers/src/KagomeSourceAudit.jl":
        "dee29180fadb1ddee104dc49a117f0447a08d100b7ebcad4efec313e45e6636b",
    "tracks/polyopt/solutions/sdp-gap-seekers/scripts/audit_exact_ray_psd.jl":
        "3fc17c8484d09bedf58d96c05a65a0f5062c804e61be00362f67ebac0f383633",
    "tracks/polyopt/solutions/sdp-gap-seekers/scripts/audit_kagome_source_assembly.jl":
        "24cd2d1b848bd1d0899c6ed4dbb50a65644cae787ca328323cb3c75f09862ac4",
    "tracks/polyopt/solutions/sdp-gap-seekers/scripts/reconstruct_exact_ray_coefficients.jl":
        "5096004b7556bd3cc82ca013be61599d775db441b911222c68ee76fa6be6955b",
    "tracks/polyopt/solutions/sdp-gap-seekers/scripts/repair_exact_psd_kernel.jl":
        "395b1e407f4cdcaf4f12bffe5165143a16035dc1dc92e738e56b2384fc4ab684",
    "tracks/polyopt/solutions/sdp-gap-seekers/spectralgap_a1171c9.patch":
        "332c0931ac810289aa3713af0948f259c01189270706af58b262d60d994d4abd",
}

SPECTRALGAP_FILES = {
    "src/SpectralGap.jl":
        "940cd72b9c4bea39b6daaefd2b9797c54df450cb0afbfb4d9322b2df1b3838bb",
    "src/basicfunction.jl":
        "2095cf7401355f37e9d17915b3ab29d44712d8e40f750eb8449f8c294229b03a",
    "src/sdp.jl":
        "4ea362723bd7601e67db3bc27f21a4a11506791ce2b9b82cb7e4a60b0f6bae10",
    "src/strengthening.jl":
        "de56b12b17049f81f689d4caef193b9dfd3bf50061fc78b1bf1547a748f7c57b",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_hash(path: Path, expected: str) -> None:
    if not path.is_file():
        raise ValueError(f"required file is absent: {path}")
    actual = sha256(path)
    if actual != expected:
        raise ValueError(
            f"SHA-256 mismatch for {path}: expected {expected}, got {actual}",
        )


def git_output(*args: str) -> str:
    return subprocess.run(
        ["git", *args],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    ).stdout.strip()


def require_line(path: Path, line: str) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    if line not in lines:
        raise ValueError(f"required audit line is absent from {path}: {line}")


def validate_exact_audit(path: Path) -> None:
    required = [
        "variable_count\t54944",
        "equality_count\t15671",
        "nonzero_exact_residuals\t0",
        "objective_improvement_exact\t119089//14841408527",
        "objective_improvement_positive\ttrue",
        "coefficient_tolerance\t1.0e-14",
        "optimizer_invoked\tfalse",
    ]
    for line in required:
        require_line(path, line)
    block_lines = [
        line for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith("block\t")
    ]
    if len(block_lines) != 9:
        raise ValueError("exact audit does not contain nine block verdicts")
    for block_index, line in enumerate(block_lines, start=1):
        if not line.startswith(f"block\t{block_index}\tproved\ttrue\t"):
            raise ValueError(f"PSD block {block_index} is not rigorously proved")
    if "exact_kernel_reduced_interval_ldlt_positive_semidefinite" not in block_lines[2]:
        raise ValueError("PSD block 3 lacks its exact-kernel proof")
    if any(
        "\tstatus\texact_zero\t" not in block_lines[index]
        for index in range(4, 9)
    ):
        raise ValueError("PSD blocks 5-9 are not all exact zero")


def validate_source_audit(path: Path) -> None:
    required = [
        f"source_commit\t{SPECTRALGAP_COMMIT}",
        "gamma\t159//125",
        "basis_dimensions\t271,104",
        "gap_basis_dimensions\t18,17",
        "strengthening_dimensions\t1,9,36,84,126",
        "stationarity_monomials\t20",
        "lambda_ordinal\t54944",
        "variable_count\t54944",
        "equality_count\t15671",
        "equality_offsets_zero\ttrue",
        "psd_dimensions\t271,104,18,17,1,9,36,84,126",
        "objective_sense\tMAX_SENSE",
        "objective\t54944 => 1//1",
        "affine_row_sha256\t7de9f007983d79adcaa0997e7a8aea170b9f73bc050d32fa6731290f11176eb2",
        "rows_compared\t15671",
        "coefficient_mismatches\t0",
        "optimizer_invoked\tfalse",
        "source_assembly_equal\ttrue",
    ]
    for line in required:
        require_line(path, line)


def relative_path(path: Path, repo_root: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(repo_root.resolve()).as_posix()
    except ValueError as error:
        raise ValueError(f"artifact is outside the repository: {path}") from error


def build_document(args: argparse.Namespace) -> dict[str, object]:
    repo_root = args.repo_root.resolve()
    git_dir = repo_root / ".bohr-handoff" / "local-git"
    if not git_dir.is_dir():
        raise ValueError("alternate local Git database is absent")
    commit = git_output(
        f"--git-dir={git_dir}",
        "rev-parse",
        f"{PROOF_COMMIT}^{{commit}}",
    )
    tree = git_output(
        f"--git-dir={git_dir}",
        "rev-parse",
        f"{PROOF_COMMIT}^{{tree}}",
    )
    if commit != PROOF_COMMIT or tree != PROOF_TREE:
        raise ValueError("proof commit/tree identity mismatch")

    for relative, expected in PROOF_FILES.items():
        require_hash(repo_root / relative, expected)

    spectralgap_root = args.spectralgap_root.resolve()
    spectralgap_commit = git_output(
        "-C",
        str(spectralgap_root),
        "rev-parse",
        "HEAD",
    )
    if spectralgap_commit != SPECTRALGAP_COMMIT:
        raise ValueError("SpectralGap commit identity mismatch")
    for relative, expected in SPECTRALGAP_FILES.items():
        require_hash(spectralgap_root / relative, expected)

    artifacts = {
        "model": args.model.resolve(),
        "ray": args.ray.resolve(),
        "exact_audit": args.exact_audit.resolve(),
        "source_audit": args.source_audit.resolve(),
    }
    for name, path in artifacts.items():
        require_hash(path, EXPECTED_ARTIFACTS[name])
    validate_exact_audit(artifacts["exact_audit"])
    validate_source_audit(artifacts["source_audit"])

    ray_line_count = sum(
        1 for _ in artifacts["ray"].open("r", encoding="utf-8")
    )
    if ray_line_count != 54_945:
        raise ValueError("exact ray does not contain 54,944 values")

    return {
        "schema": "ais-kagome-gap-certificate-v1",
        "claim": {
            "quantity": "bulk_spectral_gap",
            "relation": "less_than_or_equal",
            "gamma_numerator": 159,
            "gamma_denominator": 125,
            "gamma_decimal": "1.272",
            "scope": (
                "sign-and-cyclic-spin-symmetry-restricted "
                "infinite-volume KMS ground-state state-polynomial class"
            ),
            "unrestricted_state_claim": False,
        },
        "setup": {
            "model": "spin-half antiferromagnetic Kagome Heisenberg",
            "hamiltonian": (
                "sum_(i,j) S_i dot S_j = "
                "1/4 sum_(i,j) (X_i X_j + Y_i Y_j + Z_i Z_j)"
            ),
            "N": 13,
            "d": 3,
            "lso": 5,
            "triangles": [
                [1, 2, 3],
                [1, 4, 5],
                [2, 6, 7],
                [3, 8, 9],
                [4, 10, 11],
                [5, 12, 13],
            ],
            "inner_triangle_indices": [1, 2],
            "edges": [],
        },
        "proof": {
            "proof_commit": PROOF_COMMIT,
            "proof_tree": PROOF_TREE,
            "spectralgap_commit": SPECTRALGAP_COMMIT,
            "coefficient_tolerance": "1e-14",
            "exact_affine_rows": 15_671,
            "nonzero_exact_residuals": 0,
            "objective_improvement": "119089/14841408527",
            "psd_block_dimensions": [271, 104, 18, 17, 1, 9, 36, 84, 126],
            "psd_proof": (
                "exact zero rows/kernels plus 256-bit "
                "directed-interval LDLT"
            ),
            "source_row_sha256": (
                "7de9f007983d79adcaa0997e7a8aea170b9f73bc050d32fa6731290f11176eb2"
            ),
            "optimizer_invoked": False,
            "proof_file_sha256": PROOF_FILES,
            "spectralgap_file_sha256": SPECTRALGAP_FILES,
        },
        "artifacts": {
            name: {
                "path": relative_path(path, repo_root),
                "sha256": EXPECTED_ARTIFACTS[name],
            }
            for name, path in artifacts.items()
        },
    }


def canonical_bytes(document: dict[str, object]) -> bytes:
    return (
        json.dumps(
            document,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("utf-8")


def write_once(path: Path, payload: bytes) -> None:
    destination = path.resolve()
    if destination.exists():
        raise FileExistsError(f"refusing to overwrite certificate: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.tmp-",
        dir=destination.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--spectralgap-root", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--ray", type=Path, required=True)
    parser.add_argument("--exact-audit", type=Path, required=True)
    parser.add_argument("--source-audit", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--output", type=Path)
    mode.add_argument("--verify", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    payload = canonical_bytes(build_document(args))
    if args.output is not None:
        write_once(args.output, payload)
        print(f"certificate_path\t{args.output.resolve()}", flush=True)
        print(f"certificate_sha256\t{hashlib.sha256(payload).hexdigest()}", flush=True)
        print("certificate_materialized\ttrue", flush=True)
    else:
        actual = args.verify.read_bytes()
        if actual != payload:
            raise ValueError("certificate is not the canonical validated envelope")
        print(f"certificate_path\t{args.verify.resolve()}", flush=True)
        print(f"certificate_sha256\t{hashlib.sha256(actual).hexdigest()}", flush=True)
        print("certificate_verified\ttrue", flush=True)
    print("optimizer_invoked\tfalse", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
