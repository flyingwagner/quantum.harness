#!/usr/bin/env python3

import argparse
import json
import os
from pathlib import Path


RESULT_FIELDS = (
    "zero_threshold",
    "initial_nonzero",
    "initial_coupled_nonzero",
    "correction_nonzero",
    "correction_max_abs",
    "final_nonzero",
    "objective_improvement",
    "psd_minimum",
    "psd_minima",
    "optimizer_invoked",
)


def tolerance_key(value: str) -> str:
    return float(value).hex()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Materialize affine-projection scan cells from TSV.",
    )
    parser.add_argument("run_spec", type=Path)
    parser.add_argument("results", type=Path, nargs="+")
    return parser.parse_args()


def read_results(path: Path) -> dict[str, dict[str, object]]:
    with path.open("r", encoding="utf-8") as stream:
        header = stream.readline().rstrip("\n").split("\t")
        expected = ["rational_tolerance", *RESULT_FIELDS]
        if header != expected:
            raise ValueError("unexpected affine-projection result header")
        results: dict[str, dict[str, object]] = {}
        for line_number, line in enumerate(stream, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(expected):
                raise ValueError(f"malformed result row {line_number}")
            tolerance = tolerance_key(fields[0])
            if tolerance in results:
                raise ValueError(f"duplicate result tolerance: {tolerance}")
            result = dict(zip(RESULT_FIELDS, fields[1:], strict=True))
            result["initial_nonzero"] = int(result["initial_nonzero"])
            result["initial_coupled_nonzero"] = int(
                result["initial_coupled_nonzero"],
            )
            result["correction_nonzero"] = int(result["correction_nonzero"])
            result["final_nonzero"] = int(result["final_nonzero"])
            for key in (
                "correction_max_abs",
                "objective_improvement",
                "psd_minimum",
            ):
                result[key] = float(result[key])
            if result["optimizer_invoked"] not in ("true", "false"):
                raise ValueError("invalid optimizer_invoked value")
            result["optimizer_invoked"] = (
                result["optimizer_invoked"] == "true"
            )
            result["exact_affine_success"] = result["final_nonzero"] == 0
            results[tolerance] = result
    return results


def atomic_json(path: Path, value: object) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite manifest: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("x", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.link(temporary, path)
        temporary.unlink()
    except BaseException:
        if temporary.exists():
            temporary.unlink()
        raise


def main() -> int:
    args = parse_args()
    with args.run_spec.open("r", encoding="utf-8") as stream:
        run_spec = json.load(stream)
    results: dict[str, dict[str, object]] = {}
    for result_path in args.results:
        for tolerance, result in read_results(result_path).items():
            if tolerance in results:
                raise ValueError(
                    f"duplicate tolerance across result files: {tolerance}",
                )
            results[tolerance] = result
    expected_tolerances = {
        tolerance_key(cell["params"]["rational_tolerance"])
        for cell in run_spec["cells"]
    }
    if set(results) != expected_tolerances:
        raise ValueError("result tolerances do not match run-spec cells")
    run_directory = Path(run_spec["run_dir"])
    for cell in run_spec["cells"]:
        tolerance = tolerance_key(cell["params"]["rational_tolerance"])
        manifest = {
            "cell_id": cell["cell_id"],
            "params": cell["params"],
            "settings": {
                **run_spec["settings"],
                **cell.get("settings", {}),
            },
            "provenance": run_spec["provenance"],
            "result": results[tolerance],
        }
        atomic_json(
            run_directory / "cells" / cell["cell_id"] / "manifest.json",
            manifest,
        )
    print(f"materialized {len(run_spec['cells'])} manifests", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
