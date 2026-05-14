#!/usr/bin/env python3
"""
compare_outputs.py — KROME regression comparison.

Compares numerical files in testoutputs/ against reference files stored
inside each test's directory (tests/<test_dir>/reference/).
Writes a Markdown table to stdout and appends to $GITHUB_STEP_SUMMARY when set.
Always exits 0 (purely informational; never fails CI).

Usage
-----
    python3 compare_outputs.py                        # all tests
    python3 compare_outputs.py --test full_ismEq      # one test
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np

KROME_ROOT   = Path(__file__).resolve().parent
TESTOUTS_DIR = KROME_ROOT / "testoutputs"

# Relative deviation thresholds
THRESH_OK   = 0.01   # < 1 %   → green
THRESH_WARN = 0.10   # < 10 %  → yellow  (>= 10 % → red)

# Values smaller than this (in both ref and new) are considered negligible
# and excluded from the per-column max to avoid noise from near-zero species.
NEGLIGIBLE = 1e-30

# ---------------------------------------------------------------------------
# Test registry
#   test_dir   : sub-directory under tests/ (contains reference/ sub-dir)
#   globs      : glob patterns for numerical output files in testoutputs/
# ---------------------------------------------------------------------------
TESTS = {
    "full_ismEq": {
        "test_dir":  "popsicle_semenov_photo_cr_full_ismEqTest",
        "globs":     ["AB_Z*"],
        "skip_cols": ["t_tot", "t_cool", "n_iter"],
    },
    "full_collapse": {
        "test_dir": "popsicle_semenov_photo_cr_full_collapseTest",
        "globs":    ["fort.22"],
    },
    "full_RolligPDR": {
        "test_dir": "popsicle_semenov_photo_cr_full_RolligPDRBenchmark",
        "globs":    ["PDR_Z0"],
    },
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ref_dir(test_name: str) -> Path:
    return KROME_ROOT / "tests" / TESTS[test_name]["test_dir"] / "reference"


def _resolve_files(globs: list[str]) -> list[Path]:
    """Expand glob patterns against TESTOUTS_DIR, sorted."""
    files: list[Path] = []
    for pattern in globs:
        files.extend(sorted(TESTOUTS_DIR.glob(pattern)))
    return files


def _col_names(path: Path) -> list[str]:
    """Return column names from the first `#`-prefixed header line."""
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#"):
                return line.lstrip("#").split()
    return []


def _load(path: Path) -> np.ndarray | None:
    """Load a whitespace-separated numerical file, skipping # comments and blank lines."""
    try:
        data = np.loadtxt(path, comments="#")
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data
    except Exception:
        return None


def _max_rel_dev(new: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """
    Per-column maximum relative deviation, ignoring rows where both values
    are negligible (|new| + |ref| < NEGLIGIBLE) and NaN/Inf entries.
    """
    ncols = ref.shape[-1] if ref.ndim > 1 else 1
    if new.shape != ref.shape:
        return np.full(ncols, np.nan)

    rel = np.abs(new - ref) / (np.abs(ref) + 1e-300)

    # Mask negligible and non-finite entries
    mask = (np.abs(new) + np.abs(ref) < NEGLIGIBLE) | ~np.isfinite(rel)
    rel[mask] = np.nan

    import warnings
    with warnings.catch_warnings(), np.errstate(all="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)
        result = np.nanmax(rel, axis=0)

    return result


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------

def _emoji(dev: float) -> str:
    if np.isnan(dev):
        return "⚪"
    if dev < THRESH_OK:
        return "🟢"
    if dev < THRESH_WARN:
        return "🟡"
    return "🔴"


def _fmt_dev(dev: float) -> str:
    if np.isnan(dev):
        return "—"
    return f"{dev * 100:.2f}%"


def _compare_file(new_path: Path, ref_path: Path, skip_cols: list[str] | None = None) -> str:
    """Return a markdown section for one file comparison."""
    lines: list[str] = []
    lines.append(f"### `{new_path.name}`")

    rel_ref = ref_path.relative_to(KROME_ROOT)

    if not ref_path.exists():
        lines.append(f"> ⚪ No reference file found at `{rel_ref}` — skipping.")
        return "\n".join(lines)

    col_names = _col_names(new_path) or _col_names(ref_path)
    new = _load(new_path)
    ref = _load(ref_path)

    if new is None:
        lines.append(f"> ⚠️ Could not load `testoutputs/{new_path.name}`.")
        return "\n".join(lines)
    if ref is None:
        lines.append(f"> ⚠️ Could not load `{rel_ref}`.")
        return "\n".join(lines)

    devs = _max_rel_dev(new, ref)
    ncols = len(devs)

    # Pad / trim column names to match data width
    if len(col_names) < ncols:
        col_names += [f"col{i}" for i in range(len(col_names), ncols)]
    col_names = col_names[:ncols]

    if new.shape != ref.shape:
        lines.append(
            f"> ⚠️ Shape mismatch: new={new.shape} vs ref={ref.shape}. "
            "Per-column results shown where shapes align."
        )

    skip = set(skip_cols or [])
    lines.append("")
    lines.append("| Column | Max rel. deviation | Status |")
    lines.append("|--------|--------------------|--------|")
    for name, dev in zip(col_names, devs):
        if name not in skip:
            lines.append(f"| `{name}` | {_fmt_dev(dev)} | {_emoji(dev)} |")

    return "\n".join(lines)


def _compare_test(test_name: str) -> str:
    """Return the full markdown section for one test."""
    cfg   = TESTS[test_name]
    files = _resolve_files(cfg["globs"])
    ref_d = _ref_dir(test_name)

    sections: list[str] = [f"## Regression: {test_name}", ""]

    if not files:
        sections.append(
            f"> ⚠️ No output files found in `testoutputs/` matching `{cfg['globs']}`."
        )
        return "\n".join(sections)

    skip_cols = cfg.get("skip_cols")
    for new_path in files:
        ref_path = ref_d / new_path.name
        sections.append(_compare_file(new_path, ref_path, skip_cols=skip_cols))
        sections.append("")

    return "\n".join(sections)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1].strip())
    p.add_argument(
        "--test", choices=list(TESTS), metavar="NAME",
        help=f"compare a single test (one of: {', '.join(TESTS)})",
    )
    args = p.parse_args()

    selected = [args.test] if args.test else list(TESTS)

    output_parts: list[str] = []
    for test_name in selected:
        output_parts.append(_compare_test(test_name))

    full_output = "\n".join(output_parts)

    print(full_output)

    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary_path:
        with open(summary_path, "a") as fh:
            fh.write("\n")
            fh.write(full_output)
            fh.write("\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
