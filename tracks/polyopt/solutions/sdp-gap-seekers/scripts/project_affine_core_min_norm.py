#!/usr/bin/env python3

import argparse
import os
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Apply a floating minimum-norm projection onto an exported "
            "homogeneous affine core. This is a diagnostic, not an exact "
            "certificate."
        ),
    )
    parser.add_argument("core", type=Path)
    parser.add_argument("ray", type=Path)
    parser.add_argument("output", type=Path)
    return parser.parse_args()


def read_ray(path: Path) -> tuple[list[list[str]], np.ndarray]:
    rows: list[list[str]] = []
    values: list[float] = []
    with path.open("r", encoding="utf-8") as stream:
        if stream.readline().rstrip("\n") != "ordinal\tmoi_index\tname\tvalue":
            raise ValueError("unexpected ray header")
        for ordinal, line in enumerate(stream, start=1):
            fields = line.rstrip("\n").split("\t")
            if (
                len(fields) != 4
                or int(fields[0]) != ordinal
                or int(fields[1]) != ordinal
            ):
                raise ValueError(f"noncanonical ray row {ordinal}")
            rows.append(fields[:3])
            values.append(float(fields[3]))
    result = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(result)):
        raise ValueError("ray contains a nonfinite value")
    return rows, result


def read_core(
    path: Path,
) -> tuple[object, list[int]]:
    rows: list[int] = []
    columns: list[int] = []
    coefficients: list[float] = []
    variable_by_column: dict[int, int] = {}
    with path.open("r", encoding="utf-8") as stream:
        expected = (
            "row_position\tcolumn_position\trow_index\tvariable_index\t"
            "numerator\tdenominator"
        )
        if stream.readline().rstrip("\n") != expected:
            raise ValueError("unexpected affine-core header")
        for line_number, line in enumerate(stream, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise ValueError(f"malformed affine-core row {line_number}")
            row, column, _, variable, numerator, denominator = map(int, fields)
            if row < 1 or column < 1 or denominator <= 0:
                raise ValueError(f"invalid affine-core row {line_number}")
            previous = variable_by_column.setdefault(column, variable)
            if previous != variable:
                raise ValueError("inconsistent core variable identity")
            rows.append(row - 1)
            columns.append(column - 1)
            coefficients.append(numerator / denominator)
    row_count = max(rows) + 1
    column_count = max(columns) + 1
    if set(variable_by_column) != set(range(1, column_count + 1)):
        raise ValueError("affine core has a missing column identity")
    matrix = coo_matrix(
        (coefficients, (rows, columns)),
        shape=(row_count, column_count),
    ).tocsc()
    variables = [
        variable_by_column[column] - 1
        for column in range(1, column_count + 1)
    ]
    return matrix, variables


def write_ray(
    path: Path,
    identities: list[list[str]],
    values: np.ndarray,
) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite projected ray: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("x", encoding="utf-8") as stream:
            stream.write("ordinal\tmoi_index\tname\tvalue\n")
            for identity, value in zip(identities, values, strict=True):
                stream.write("\t".join([*identity, repr(float(value))]) + "\n")
        os.link(temporary, path)
        temporary.unlink()
    except BaseException:
        if temporary.exists():
            temporary.unlink()
        raise


def main() -> int:
    args = parse_args()
    identities, ray = read_ray(args.ray)
    core, variables = read_core(args.core)
    scale = np.max(np.abs(ray))
    if not scale > 0:
        raise ValueError("ray is identically zero")
    normalized = ray / scale
    core_values = normalized[variables]
    initial_residual = core @ core_values
    normal = (core @ core.T).tocsc()
    factor = splu(normal, permc_spec="COLAMD")
    multipliers = factor.solve(-initial_residual)
    correction = core.T @ multipliers
    projected_core = core_values + correction
    final_residual = core @ projected_core
    projected = normalized.copy()
    projected[variables] = projected_core
    write_ray(args.output, identities, projected)

    print(f"rows\t{core.shape[0]}", flush=True)
    print(f"columns\t{core.shape[1]}", flush=True)
    print(f"normal_nonzeros\t{normal.nnz}", flush=True)
    print(f"factor_l_nonzeros\t{factor.L.nnz}", flush=True)
    print(f"factor_u_nonzeros\t{factor.U.nnz}", flush=True)
    print(
        f"initial_residual_infinity\t"
        f"{np.linalg.norm(initial_residual, np.inf)}",
        flush=True,
    )
    print(
        f"final_residual_infinity\t"
        f"{np.linalg.norm(final_residual, np.inf)}",
        flush=True,
    )
    print(
        f"correction_infinity\t{np.linalg.norm(correction, np.inf)}",
        flush=True,
    )
    print("exact_projection\tfalse", flush=True)
    print("optimizer_invoked\tfalse", flush=True)
    print(f"output\t{args.output.resolve()}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
