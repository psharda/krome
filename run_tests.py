#!/usr/bin/env python3
"""
run_tests.py — KROME chemical network test suite runner.

Runs code generation, compilation, execution, and plotting for the tests
described in testSuite.MD.  Timing is recorded for each phase.  After all
tests complete, output files are collected into Outputs/<build_name>/.

Usage
-----
    python3 run_tests.py                                 # run all tests
    python3 run_tests.py --tests full_ismEq gow_ismEq   # named subset
    python3 run_tests.py --list                          # list test names
    python3 run_tests.py --skip-generate                 # reuse existing builds
    python3 run_tests.py --skip-compile                  # skip compilation
    python3 run_tests.py --skip-run                      # skip running
    python3 run_tests.py --skip-plot                     # skip plotting
    python3 run_tests.py --no-collect                    # skip output collection
"""

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

KROME_ROOT = Path(__file__).resolve().parent
OUTPUTS_DIR = KROME_ROOT / "testoutputs"

# ---------------------------------------------------------------------------
# Test registry
# ---------------------------------------------------------------------------
# Each entry describes one KROME test:
#   test_dir     : sub-directory under tests/
#   build_name   : label passed to -project; build dir = build_<build_name>/
#   executables  : list of compiled binaries to run, in order
#   output_globs : glob patterns (relative to build dir) for data output files
#   plot_script  : Python plotting script to run after executables; None = skip
#   plot_outputs : glob patterns for files created by the plot script
# ---------------------------------------------------------------------------
TESTS = [
    {
        "test_dir":    "popsicle_semenov_photo_cr_full_ismEqTest",
        "build_name":  "full_ismEq",
        "executables": ["test_cooling", "test_shielded"],
        "output_globs": ["AB_Z*"],
        "plot_script": "plot.py",
        "plot_outputs": ["ISMEq_Z.pdf"],
    },
    {
        "test_dir":    "popsicle_semenov_photo_cr_full_collapseTest",
        "build_name":  "full_collapse",
        "executables": ["test"],
        "output_globs": ["fort.22", "fort.31", "fort.911"],
        "plot_script": "plot.py",
        "plot_outputs": ["CollapseTgas_Metallicities*"],
    },
    {
        "test_dir":    "popsicle_semenov_photo_cr_full_RolligPDRBenchmark",
        "build_name":  "full_RolligPDR",
        "executables": ["test_pdr"],
        "output_globs": ["PDR_Z0"],
        "plot_script": "plot.py",
        "plot_outputs": ["RolligPDR.pdf"],
    },
]

TEST_BY_NAME = {t["build_name"]: t for t in TESTS}

# ---------------------------------------------------------------------------
# Terminal helpers
# ---------------------------------------------------------------------------

RESET  = "\033[0m"
BOLD   = "\033[1m"
GREEN  = "\033[32m"
YELLOW = "\033[33m"
RED    = "\033[31m"
CYAN   = "\033[36m"


def _header(msg: str) -> None:
    bar = "=" * 66
    print(f"\n{BOLD}{CYAN}{bar}\n  {msg}\n{bar}{RESET}")


def _step(msg: str) -> None:
    print(f"\n{BOLD}▶  {msg}{RESET}")


def _ok(msg: str) -> None:
    print(f"{GREEN}✓  {msg}{RESET}")


def _warn(msg: str) -> None:
    print(f"{YELLOW}⚠  {msg}{RESET}", file=sys.stderr)


def _err(msg: str) -> None:
    print(f"{RED}✗  {msg}{RESET}", file=sys.stderr)


def _fmt(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    return f"{int(seconds // 60)}m {seconds % 60:.0f}s"


# ---------------------------------------------------------------------------
# Command runner
# ---------------------------------------------------------------------------

def _run(cmd: list, cwd: Path) -> tuple[int, float]:
    """Run *cmd* in *cwd*, streaming stdout/stderr. Returns (rc, elapsed_s)."""
    t0 = time.perf_counter()
    rc = subprocess.run(cmd, cwd=cwd).returncode
    return rc, time.perf_counter() - t0


# ---------------------------------------------------------------------------
# Phases
# ---------------------------------------------------------------------------

def generate(test: dict) -> tuple[bool, float]:
    """Invoke krome to produce Fortran source in build_<build_name>/.

    Newlines are piped to stdin so that krome's interactive warning prompts
    (e.g. "Any key to continue q to quit...") are dismissed automatically.
    """
    _step(f"[generate]  krome -test={test['test_dir']} -project={test['build_name']}")
    cmd = [sys.executable, "krome",
           f"-test={test['test_dir']}",
           f"-project={test['build_name']}"]
    t0 = time.perf_counter()
    # Pipe enough newlines to answer any number of interactive warning prompts.
    result = subprocess.run(cmd, cwd=KROME_ROOT, input=b"\n" * 20)
    t = time.perf_counter() - t0
    if result.returncode != 0:
        _err(f"Code generation failed (exit {result.returncode})")
        return False, t
    _ok(f"Generated  ({_fmt(t)})")
    return True, t


def compile_test(test: dict) -> tuple[bool, float]:
    """Compile the generated code with gfortran."""
    build_dir = KROME_ROOT / f"build_{test['build_name']}"
    if not build_dir.is_dir():
        _err(f"Build directory not found: {build_dir.name}/  (run without --skip-generate)")
        return False, 0.0
    _step(f"[compile]   make gfortran  ({build_dir.name}/)")
    rc, t = _run(["make", "gfortran"], cwd=build_dir)
    if rc != 0:
        _err(f"Compilation failed (exit {rc})")
        return False, t
    _ok(f"Compiled  ({_fmt(t)})")
    return True, t


def run_executables(test: dict) -> tuple[bool, dict]:
    """Run each listed executable in order; return (all_passed, {exe: elapsed})."""
    build_dir = KROME_ROOT / f"build_{test['build_name']}"
    timings: dict[str, float] = {}
    all_passed = True

    for exe in test["executables"]:
        exe_path = build_dir / exe
        if not exe_path.exists():
            _err(f"Executable not found: {exe_path}  (compilation may have failed)")
            timings[exe] = 0.0
            all_passed = False
            continue
        _step(f"[run]       ./{exe}  ({build_dir.name}/)")
        rc, t = _run([f"./{exe}"], cwd=build_dir)
        timings[exe] = t
        if rc != 0:
            _err(f"./{exe} failed (exit {rc}) after {_fmt(t)}")
            all_passed = False
        else:
            _ok(f"./{exe}  ({_fmt(t)})")

    return all_passed, timings


def run_plot(test: dict) -> tuple[bool, float]:
    """Run the plot script from inside the build directory."""
    script = test.get("plot_script")
    if not script:
        return True, 0.0
    build_dir = KROME_ROOT / f"build_{test['build_name']}"
    if not (build_dir / script).exists():
        _warn(f"Plot script {script} not found in {build_dir.name}/; skipping")
        return True, 0.0
    _step(f"[plot]      python3 {script}  ({build_dir.name}/)")
    rc, t = _run([sys.executable, script], cwd=build_dir)
    if rc != 0:
        _warn(f"Plot script exited {rc} — outputs may be incomplete  ({_fmt(t)})")
        return False, t
    _ok(f"Plots written  ({_fmt(t)})")
    return True, t


# ---------------------------------------------------------------------------
# Output collection
# ---------------------------------------------------------------------------

def collect(test: dict) -> None:
    """Move plot-input data and PDFs from build_<name>/ into testoutputs/."""
    build_dir = KROME_ROOT / f"build_{test['build_name']}"
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)

    patterns = test.get("output_globs", []) + test.get("plot_outputs", [])
    moved: list[str] = []

    for pattern in patterns:
        for src in sorted(build_dir.glob(pattern)):
            dest = OUTPUTS_DIR / src.name
            shutil.move(str(src), dest)
            moved.append(src.name)

    if moved:
        _ok(f"Moved {len(moved)} file(s)  →  testoutputs/")
        for name in moved:
            print(f"      {name}")
    else:
        _warn(f"No output files matched for {test['build_name']}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__.split("\n")[1].strip(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Available tests:\n" + "\n".join(
            f"  {t['build_name']:<28} ({t['test_dir']})" for t in TESTS
        ),
    )
    p.add_argument(
        "--tests", nargs="+", metavar="NAME",
        help="run only the named test(s); default: all (see --list)",
    )
    p.add_argument("--list", action="store_true",
                   help="print available test names and exit")
    p.add_argument("--skip-generate", action="store_true",
                   help="skip code generation (reuse existing build directory)")
    p.add_argument("--skip-compile",  action="store_true",
                   help="skip compilation")
    p.add_argument("--skip-run",      action="store_true",
                   help="skip running executables")
    p.add_argument("--skip-plot",     action="store_true",
                   help="skip plotting")
    p.add_argument("--no-collect",    action="store_true",
                   help="skip moving outputs to testoutputs/")
    return p.parse_args()


def main() -> int:
    args = _parse_args()

    if args.list:
        print("Available tests (pass one or more names to --tests):")
        for t in TESTS:
            print(f"  {t['build_name']:<28} ({t['test_dir']})")
        return 0

    # resolve selected tests
    if args.tests:
        unknown = [n for n in args.tests if n not in TEST_BY_NAME]
        if unknown:
            _err(f"Unknown test name(s): {', '.join(unknown)}")
            _err(f"Run with --list to see available names.")
            return 1
        selected = [TEST_BY_NAME[n] for n in args.tests]
    else:
        selected = TESTS

    # timing summary accumulated per test
    # { build_name: {"generate": s, "compile": s, "run": {exe: s}, "plot": s} }
    summary: dict[str, dict] = {}
    failed:  list[str] = []

    for test in selected:
        build_name = test["build_name"]
        _header(f"TEST: {build_name}")
        summary[build_name] = {"generate": 0.0, "compile": 0.0, "run": {}, "plot": 0.0}
        test_passed = True

        if not args.skip_generate:
            ok_gen, t = generate(test)
            summary[build_name]["generate"] = t
            if not ok_gen:
                failed.append(build_name)
                continue   # cannot compile or run without generated source

        if not args.skip_compile:
            ok_cmp, t = compile_test(test)
            summary[build_name]["compile"] = t
            if not ok_cmp:
                failed.append(build_name)
                continue   # cannot run without compiled binaries

        if not args.skip_run:
            ok_run, run_times = run_executables(test)
            summary[build_name]["run"] = run_times
            if not ok_run:
                test_passed = False

        if not args.skip_plot:
            ok_plt, t = run_plot(test)
            summary[build_name]["plot"] = t
            if not ok_plt:
                test_passed = False

        if not test_passed:
            failed.append(build_name)

    # collect outputs after all tests have run
    if not args.no_collect:
        _header("Collecting outputs  →  testoutputs/")
        for test in selected:
            collect(test)

    # timing summary table
    _header("Timing summary")
    col = 26
    print(f"{'Test':<{col}}  {'Generate':>9}  {'Compile':>9}  {'Run':>12}  {'Plot':>8}  {'Total':>9}")
    print("-" * 80)
    for build_name, times in summary.items():
        gen_t  = times.get("generate", 0.0)
        cmp_t  = times.get("compile",  0.0)
        run_t  = sum(times.get("run", {}).values())
        plt_t  = times.get("plot",     0.0)
        total  = gen_t + cmp_t + run_t + plt_t
        status = f"{RED}FAILED{RESET}" if build_name in failed else f"{GREEN}ok{RESET}"
        print(
            f"{build_name:<{col}}  "
            f"{_fmt(gen_t):>9}  {_fmt(cmp_t):>9}  "
            f"{_fmt(run_t):>12}  {_fmt(plt_t):>8}  "
            f"{_fmt(total):>9}  {status}"
        )

    if failed:
        print(f"\n{RED}Failed: {', '.join(failed)}{RESET}")
        return 1

    print(f"\n{GREEN}All {len(selected)} test(s) completed successfully.{RESET}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
