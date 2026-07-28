#!/usr/bin/env python3

import argparse
import gc
import os
from pathlib import Path

from flint import fmpz, nmod_mat


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rank an exact rational affine minor over prime fields.",
    )
    parser.add_argument("minor", type=Path)
    parser.add_argument("primes", type=int, nargs="+")
    parser.add_argument(
        "--pivots-output",
        type=Path,
        help="write pivot columns from RREF over the first prime",
    )
    return parser.parse_args()


def read_matrix(
    path: Path,
) -> tuple[
    int,
    int,
    list[tuple[int, int, int, int]],
    list[int],
]:
    entries: list[tuple[int, int, int, int]] = []
    variable_by_column: dict[int, int] = {}
    row_count = 0
    column_count = 0
    with path.open("r", encoding="utf-8") as stream:
        header = stream.readline().rstrip("\n")
        expected = (
            "row_position\tcolumn_position\trow_index\tvariable_index\t"
            "numerator\tdenominator"
        )
        if header != expected:
            raise ValueError("unexpected affine-minor header")
        for line_number, line in enumerate(stream, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise ValueError(f"malformed affine-minor row {line_number}")
            row, column, _, variable, numerator, denominator = map(
                int,
                fields,
            )
            if row < 1 or column < 1 or denominator <= 0:
                raise ValueError(f"invalid affine-minor row {line_number}")
            entries.append((row - 1, column - 1, numerator, denominator))
            previous = variable_by_column.setdefault(column, variable)
            if previous != variable:
                raise ValueError(
                    f"inconsistent variable identity in column {column}",
                )
            row_count = max(row_count, row)
            column_count = max(column_count, column)
    if row_count == 0 or column_count == 0:
        raise ValueError("affine minor is empty")
    if set(variable_by_column) != set(range(1, column_count + 1)):
        raise ValueError("affine matrix has a missing column identity")
    return (
        row_count,
        column_count,
        entries,
        [variable_by_column[column] for column in range(1, column_count + 1)],
    )


def modular_rank(
    row_count: int,
    column_count: int,
    entries: list[tuple[int, int, int, int]],
    prime: int,
) -> int:
    if not fmpz(prime).is_prime():
        raise ValueError(f"modulus is not prime: {prime}")
    matrix = nmod_mat(row_count, column_count, prime)
    for row, column, numerator, denominator in entries:
        denominator_mod = denominator % prime
        if denominator_mod == 0:
            raise ValueError(
                f"coefficient denominator is zero modulo {prime}",
            )
        matrix[row, column] = (
            (numerator % prime) * pow(denominator_mod, -1, prime)
        ) % prime
    rank = matrix.rank()
    del matrix
    gc.collect()
    return rank


def modular_rref_pivots(
    row_count: int,
    column_count: int,
    entries: list[tuple[int, int, int, int]],
    prime: int,
) -> tuple[int, list[int]]:
    if not fmpz(prime).is_prime():
        raise ValueError(f"modulus is not prime: {prime}")
    matrix = nmod_mat(row_count, column_count, prime)
    for row, column, numerator, denominator in entries:
        denominator_mod = denominator % prime
        if denominator_mod == 0:
            raise ValueError(
                f"coefficient denominator is zero modulo {prime}",
            )
        matrix[row, column] = (
            (numerator % prime) * pow(denominator_mod, -1, prime)
        ) % prime
    _, rank = matrix.rref(inplace=True)
    pivots: list[int] = []
    search_start = 0
    for row in range(rank):
        pivot = next(
            (
                column
                for column in range(search_start, column_count)
                if matrix[row, column] != 0
            ),
            None,
        )
        if pivot is None:
            raise RuntimeError(f"RREF row {row + 1} has no pivot")
        pivots.append(pivot)
        search_start = pivot + 1
    del matrix
    gc.collect()
    return rank, pivots


def write_pivots(
    path: Path,
    pivots: list[int],
    variable_by_column: list[int],
    prime: int,
) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite pivots: {path}")
    temporary = path.with_name(f".{path.name}.tmp")
    if temporary.exists():
        raise FileExistsError(f"temporary pivot path exists: {temporary}")
    try:
        with temporary.open("x", encoding="utf-8") as stream:
            stream.write(
                "pivot_ordinal\tcolumn_position\tvariable_index\tprime\n",
            )
            for ordinal, column in enumerate(pivots, start=1):
                stream.write(
                    f"{ordinal}\t{column + 1}\t"
                    f"{variable_by_column[column]}\t{prime}\n",
                )
        os.link(temporary, path)
        temporary.unlink()
    except BaseException:
        if temporary.exists():
            temporary.unlink()
        raise


def main() -> int:
    args = parse_args()
    row_count, column_count, entries, variable_by_column = read_matrix(
        args.minor,
    )
    print(f"rows\t{row_count}", flush=True)
    print(f"columns\t{column_count}", flush=True)
    print(f"nonzero_coefficients\t{len(entries)}", flush=True)
    all_full = True
    for prime_index, prime in enumerate(args.primes):
        print(f"progress\tprime\t{prime}", flush=True)
        if prime_index == 0 and args.pivots_output is not None:
            rank, pivots = modular_rref_pivots(
                row_count,
                column_count,
                entries,
                prime,
            )
            if rank == row_count:
                write_pivots(
                    args.pivots_output,
                    pivots,
                    variable_by_column,
                    prime,
                )
                print(
                    f"pivots_output\t{args.pivots_output.resolve()}",
                    flush=True,
                )
            else:
                print("pivots_output\tunavailable_rank_deficient", flush=True)
        else:
            rank = modular_rank(
                row_count,
                column_count,
                entries,
                prime,
            )
        print(f"rank_mod_{prime}\t{rank}", flush=True)
        full = rank == row_count
        print(
            f"full_row_rank_mod_{prime}\t{str(full).lower()}",
            flush=True,
        )
        all_full = all_full and full
    print(f"all_primes_full_row_rank\t{str(all_full).lower()}", flush=True)
    return 0 if all_full else 2


if __name__ == "__main__":
    raise SystemExit(main())
