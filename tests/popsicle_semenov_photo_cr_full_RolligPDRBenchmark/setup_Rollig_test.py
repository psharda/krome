#!/usr/bin/env python3
"""
Update krome_ode.f90 so that the H2 formation on dust grains uses the
expression expected by the Rollig PDR benchmark tests by default.

Specifically, replace any line containing the substring:
  "nH2dust = nH2dust + "
with exactly:
  "nH2dust = nH2dust + 3d-18*sqrt(Tgas)*n(idx_H)*nH*dust2gas_ratio"

Usage:
  - No args: patches ./krome_ode.f90 in-place (default for Makefile)
  - With a path: python3 scripts/update_krome_ode.py path/to/file.f90

This default is necessary because the Rollig tests assume this rate
for H2 formation on grains.
"""

from __future__ import annotations

import sys
from pathlib import Path


TARGET_SUBSTR = "nH2dust = nH2dust + "
REPLACEMENT = (
    "nH2dust = nH2dust + 3d-18*sqrt(Tgas)*n(idx_H)*nH*dust2gas_ratio"
)


def replace_lines(file_path: Path) -> int:
    text = file_path.read_text(encoding="utf-8")
    lines = text.splitlines(True)  # keep line endings
    changed = 0
    new_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]
        # Check multi-line case first
        if (
            TARGET_SUBSTR in line.replace("&", "")
            and line.rstrip().endswith("&")
            and i + 1 < len(lines)
        ):
            indent = line[: len(line) - len(line.lstrip(" \t"))]
            newline = "\r\n" if line.endswith("\r\n") else "\n"
            new_lines.append(f"{indent}{REPLACEMENT}{newline}")
            changed += 1
            i += 2  # Skip both lines
        elif TARGET_SUBSTR in line:
            indent = line[: len(line) - len(line.lstrip(" \t"))]
            newline = "\r\n" if line.endswith("\r\n") else "\n"
            new_lines.append(f"{indent}{REPLACEMENT}{newline}")
            changed += 1
            i += 1
        else:
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

