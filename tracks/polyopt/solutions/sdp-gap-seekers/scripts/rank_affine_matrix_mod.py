#!/usr/bin/env python3

import argparse
import gc
from pathlib import Path

from flint import fmpz, nmod_mat


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rank an exact rational affine minor over prime fields.",
    )
    parser.add_argument("minor", type=Path)
    parser.add_argument("primes", type=int, nargs="+")
    return parser.parse_args()


def read_matrix(
    path: Path,
) -> tuple[int, int, list[tuple[int, int, int, int]]]:
    entries: list[tuple[int, int, int, int]] = []
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
            row, column, _, _, numerator, denominator = map(int, fields)
            if row < 1 or column < 1 or denominator <= 0:
                raise ValueError(f"invalid affine-minor row {line_number}")
            entries.append((row - 1, column - 1, numerator, denominator))
            row_count = max(row_count, row)
            column_count = max(column_count, column)
    if row_count == 0 or column_count == 0:
        raise ValueError("affine minor is empty")
    return row_count, column_count, entries


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


def main() -> int:
    args = parse_args()
    row_count, column_count, entries = read_matrix(args.minor)
    print(f"rows\t{row_count}", flush=True)
    print(f"columns\t{column_count}", flush=True)
    print(f"nonzero_coefficients\t{len(entries)}", flush=True)
    all_full = True
    for prime in args.primes:
        print(f"progress\tprime\t{prime}", flush=True)
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
