#!/usr/bin/env python3
"""
Update krome_ode.f90 so that:

1. H2 formation on dust grains uses the expression expected by the
   Rollig PDR benchmark tests by default.

   Replace any line containing:
       "nH2dust = nH2dust + "
   with:
       "nH2dust = nH2dust + 3d-18*sqrt(Tgas)*n(idx_H)*nH*dust2gas_ratio"

2. HD formation on dust grains is disabled.

   Replace any line containing:
       "nHDdust = nHDdust + "
   with:
       "nHDdust = 0d0"

Usage:

- No args: patches ./krome_ode.f90 in-place (default for Makefile)
- With a path:
      python3 scripts/update_krome_ode.py path/to/file.f90
"""

from __future__ import annotations

import sys
from pathlib import Path

REPLACEMENTS = {
    "nH2dust = nH2dust + ":
        "nH2dust = nH2dust + 3d-18*sqrt(Tgas)*n(idx_H)*nH*dust2gas_ratio",

    "nHDdust = nHDdust + ":
        "nHDdust = 0d0",
}


def replace_lines(file_path: Path) -> int:
    text = file_path.read_text(encoding="utf-8")
    lines = text.splitlines(True)  # keep line endings

    changed = 0
    new_lines = []
    i = 0

    while i < len(lines):
        line = lines[i]
        matched = False

        for target_substr, replacement in REPLACEMENTS.items():

            # Handle two-line continuations:
            #   nH2dust = nH2dust + ... &
            #   ...
            if (
                target_substr in line.replace("&", "")
                and line.rstrip().endswith("&")
                and i + 1 < len(lines)
            ):
                indent = line[: len(line) - len(line.lstrip(" \t"))]
                newline = "\r\n" if line.endswith("\r\n") else "\n"

                new_lines.append(f"{indent}{replacement}{newline}")
                changed += 1
                i += 2
                matched = True
                break

            # Handle single-line statements
            if target_substr in line:
                indent = line[: len(line) - len(line.lstrip(" \t"))]
                newline = "\r\n" if line.endswith("\r\n") else "\n"

                new_lines.append(f"{indent}{replacement}{newline}")
                changed += 1
                i += 1
                matched = True
                break

        if not matched:
            new_lines.append(line)
            i += 1

    if changed:
        file_path.write_text("".join(new_lines), encoding="utf-8")

    return changed


def main(argv: list[str]) -> int:
    # Default to local krome_ode.f90 to align with Makefile
    target = Path("krome_ode.f90") if len(argv) == 1 else Path(argv[1])

    if not target.is_file():
        print(f"Error: file not found: {target}")
        return 1

    count = replace_lines(target)
    print(f"Replacements applied: {count}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
