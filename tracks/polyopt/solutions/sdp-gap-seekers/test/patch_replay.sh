#!/usr/bin/env bash
set -euo pipefail

BASE_COMMIT="a1171c906ff2cc2901e58c2426397a2f68c32bb7"
BASE_TREE="52d2b037c2275866d482a7dc531198d412c566e1"
PATCH_SHA256="5ef9585c71b84b7a07b36610e2bc8aab060a40a8b5062633b070c92dc74fc947"
SPECTRALGAP_SHA256="940cd72b9c4bea39b6daaefd2b9797c54df450cb0afbfb4d9322b2df1b3838bb"
BASICFUNCTION_SHA256="2095cf7401355f37e9d17915b3ab29d44712d8e40f750eb8449f8c294229b03a"
SDP_SHA256="b1fa2280cca51fca38154daf5c767f7538ab68c2297e673eef474da3505f0ccc"
STRENGTHENING_SHA256="de56b12b17049f81f689d4caef193b9dfd3bf50061fc78b1bf1547a748f7c57b"

TEST_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
SOLUTION_DIR="$(dirname -- "$TEST_DIR")"
PATCH_PATH="$SOLUTION_DIR/spectralgap_a1171c9.patch"
SOURCE_REPO="${1:-https://github.com/wangjie212/SpectralGap.git}"
WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/spectralgap-patch-replay.XXXXXX")"
trap 'rm -rf -- "$WORK_DIR"' EXIT

[[ "$(sha256sum "$PATCH_PATH" | awk '{print $1}')" == "$PATCH_SHA256" ]]
git -C "$SOURCE_REPO" cat-file -e "$BASE_COMMIT^{commit}" 2>/dev/null ||
  git clone --quiet --no-checkout "$SOURCE_REPO" "$WORK_DIR/upstream"
if [[ -d "$WORK_DIR/upstream/.git" ]]; then
  SOURCE_REPO="$WORK_DIR/upstream"
fi
[[ "$(git -C "$SOURCE_REPO" rev-parse "$BASE_COMMIT^{tree}")" == "$BASE_TREE" ]]

git clone --quiet --no-checkout "$SOURCE_REPO" "$WORK_DIR/replay"
git -C "$WORK_DIR/replay" checkout --quiet --detach "$BASE_COMMIT"
git -C "$WORK_DIR/replay" apply --check "$PATCH_PATH"
git -C "$WORK_DIR/replay" apply "$PATCH_PATH"
git -C "$WORK_DIR/replay" diff --check

check_sha256() {
  local expected="$1"
  local path="$2"
  [[ "$(sha256sum "$path" | awk '{print $1}')" == "$expected" ]]
}

check_sha256 "$SPECTRALGAP_SHA256" "$WORK_DIR/replay/src/SpectralGap.jl"
check_sha256 "$BASICFUNCTION_SHA256" "$WORK_DIR/replay/src/basicfunction.jl"
check_sha256 "$SDP_SHA256" "$WORK_DIR/replay/src/sdp.jl"
check_sha256 "$STRENGTHENING_SHA256" "$WORK_DIR/replay/src/strengthening.jl"

python3 - "$WORK_DIR/replay" <<'PY'
import pathlib
import re
import sys

root = pathlib.Path(sys.argv[1])
project_text = (root / "Project.toml").read_text(encoding="utf-8")
module_text = (root / "src" / "SpectralGap.jl").read_text(encoding="utf-8")
deps_match = re.search(
    r"^\[deps\]\s*$\n(?P<body>.*?)(?=^\[|\Z)",
    project_text,
    re.MULTILINE | re.DOTALL,
)
if deps_match is None:
    raise SystemExit("SpectralGap Project.toml has no [deps] section")
declared = {
    match.group(1)
    for match in re.finditer(
        r"^([A-Za-z_][A-Za-z0-9_]*)\s*=",
        deps_match.group("body"),
        re.MULTILINE,
    )
}
imports = {
    match.group(1)
    for match in re.finditer(
        r"^\s*(?:using|import)\s+([A-Za-z_][A-Za-z0-9_]*)",
        module_text,
        re.MULTILINE,
    )
}
undeclared = imports.difference(declared, {"Base", "Core"})
if undeclared:
    raise SystemExit(f"undeclared SpectralGap imports: {sorted(undeclared)}")
if "Clarabel" in imports or "_select_optimizer" in module_text:
    raise SystemExit("obsolete Clarabel selector survived patch replay")
PY

if git -C "$WORK_DIR/replay" diff --unified=0 |
  grep -Eq '^\+\s*(using|import)\s'; then
  echo "patch unexpectedly adds a package import" >&2
  exit 1
fi

printf 'PASS patch=%s base=%s tree=%s source_hashes=4 dependency_coherence=ok\n' \
  "$PATCH_SHA256" "$BASE_COMMIT" "$BASE_TREE"
