#!/usr/bin/env python3
"""
Run the 4-version Rollig PDR benchmark suite (Rollig et al. 2007, A&A, 467, 187)
and produce a combined comparison plot against the published reference models.

Usage (from the build directory):
    python3 run_Rollig_suite.py

Workflow:
  1. For each version, patch test_pdr.f90 with the correct ntot/chi0/output
     filenames, compile a named binary (test_pdr_V1 ... test_pdr_V4), then
     restore the original source.  Compilations are sequential because they
     share a single source file.
  2. All four binaries are launched in parallel; stdout/stderr go to
     run_V1.log ... run_V4.log.
  3. Once all runs finish, plot_figure4() overlays KROME results on the
     published Rollig reference curves and saves RolligCombined.pdf.

Benchmark versions:
  V1: ntot = 10^3  cm^-3, chi = 10
  V2: ntot = 10^3  cm^-3, chi = 10^5
  V3: ntot = 10^5.5 cm^-3, chi = 10
  V4: ntot = 10^5.5 cm^-3, chi = 10^5
"""

import subprocess
import shutil
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import glob

BUILD_DIR = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(BUILD_DIR, "test_pdr.f90")

# (version label, chi0 Fortran literal, log10(ntot) Fortran literal)
VERSIONS = [
    ("V1", "1d1",  "3d0"),
    ("V2", "1d5",  "3d0"),
    ("V3", "1d1",  "5.5d0"),
    ("V4", "1d5",  "5.5d0"),
]

# Exact Fortran block that derives output filenames from metallicity in the
# original test_pdr.f90.  We replace it per-version with hardcoded names
# (PDR_V1, COOL_V1, HEAT_V1, …) so all four binaries can run simultaneously
# without clobbering each other's output files.
ORIGINAL_FILENAME_BLOCK = """\
    !Deduce filename from metallicity
    zint = int(log10(zs(jz2)))
    if(zint == 0) then
      write(zint_str, '(I1)') zint
    else
      write(zint_str, '(I2)') zint
    endif
    filename = trim('PDR_Z') // trim(zint_str)
    filename = trim(filename)
    !Open file
    open(unit=22,file=filename,status='replace',action='write')
    write(22, '(A)', ADVANCE='NO') "#ntot rho Tgas Tdust ColumnTot"
    write(22, '(A)', ADVANCE='NO') trim(krome_get_names_header())
    write(22, '(A)') " t_tot t_cool n_iter"

    filename = trim('COOL_Z') // trim(zint_str)
    filename = trim(filename)
    open(unit=31,file=filename,status='replace',action='write')
    write(31, '(A)', ADVANCE='NO') "#ColumnTot Tgas sum(cools)"
    write(31, '(A)') trim(krome_get_cooling_names_header())

    filename = trim('HEAT_Z') // trim(zint_str)
    filename = trim(filename)
    open(unit=911,file=filename,status='replace',action='write')
    write(911, '(A)', ADVANCE='NO') "#ColumnTot Tgas sum(heats)"
    write(911, '(A)') trim(krome_get_heating_names_header())"""


def make_filename_block(ver):
    return f"""\
    !Hardcoded filename for benchmark version {ver}
    open(unit=22,file='PDR_{ver}',status='replace',action='write')
    write(22, '(A)', ADVANCE='NO') "#ntot rho Tgas Tdust ColumnTot"
    write(22, '(A)', ADVANCE='NO') trim(krome_get_names_header())
    write(22, '(A)') " t_tot t_cool n_iter"

    open(unit=31,file='COOL_{ver}',status='replace',action='write')
    write(31, '(A)', ADVANCE='NO') "#ColumnTot Tgas sum(cools)"
    write(31, '(A)') trim(krome_get_cooling_names_header())

    open(unit=911,file='HEAT_{ver}',status='replace',action='write')
    write(911, '(A)', ADVANCE='NO') "#ColumnTot Tgas sum(heats)"
    write(911, '(A)') trim(krome_get_heating_names_header())"""


def patch_source(ver, chi0_val, ntot_val):
    with open(SRC) as f:
        original_text = f.read()
    text = original_text.replace('\r\n', '\n')

    # A: chi0
    text = text.replace("  chi0 = 1d1\n", f"  chi0 = {chi0_val}\n", 1)
    # B: ntot
    text = text.replace(
        "    ntot = 10**(3d0)    ! Fixed density of 100cm^-3\n",
        f"    ntot = 10**({ntot_val})  ! {ver} density\n", 1,
    )
    # C: filename block
    text = text.replace(ORIGINAL_FILENAME_BLOCK, make_filename_block(ver), 1)

    normalized_original = original_text.replace('\r\n', '\n')
    assert text != normalized_original, \
        f"[{ver}] No changes applied — check literal strings in patch_source"

    with open(SRC, 'w') as f:
        f.write(text)
    return original_text


def compile_version(ver, chi0_val, ntot_val):
    original_text = patch_source(ver, chi0_val, ntot_val)
    # Force recompilation of test_pdr.f90 by removing the stale object file.
    # Without this, make can skip the compile step when all four patches happen
    # within the same filesystem timestamp granularity (1 s on HFS+).
    obj = os.path.join(BUILD_DIR, "test_pdr.o")
    if os.path.exists(obj):
        os.remove(obj)

    print(f"[{ver}] compiling (chi0={chi0_val}, ntot=10**{ntot_val})...")
    try:
        result = subprocess.run(
            ["make", "gfortran"],
            cwd=BUILD_DIR,
            capture_output=True, text=True, check=True,
        )
        print(f"[{ver}] compile OK")
    except subprocess.CalledProcessError as e:
        print(f"[{ver}] COMPILE FAILED:\n{e.stdout}\n{e.stderr}")
        raise
    finally:
        # Always restore source, even on compile failure, so the repo stays clean.
        with open(SRC, 'w') as f:
            f.write(original_text)

    src_bin = os.path.join(BUILD_DIR, "test_pdr")
    dest_bin = os.path.join(BUILD_DIR, f"test_pdr_{ver}")
    shutil.copy2(src_bin, dest_bin)
    os.chmod(dest_bin, 0o755)
    print(f"[{ver}] binary → test_pdr_{ver}")


def run_all_parallel():
    # Each binary writes to its own output files (PDR_V1, COOL_V1, …), so
    # running all four from the same working directory is safe.
    procs = {}
    logs = {}
    for ver, _, _ in VERSIONS:
        logpath = os.path.join(BUILD_DIR, f"run_{ver}.log")
        logfile = open(logpath, 'w')
        logs[ver] = logfile
        proc = subprocess.Popen(
            [f"./test_pdr_{ver}"],
            cwd=BUILD_DIR,          # data files (.dat, .gfe) are resolved here
            stdout=logfile,
            stderr=subprocess.STDOUT,
        )
        procs[ver] = proc
        print(f"[{ver}] started PID {proc.pid}, log: run_{ver}.log")
    return procs, logs


def wait_and_check(procs, logs):
    failed = []
    for ver, proc in procs.items():
        proc.wait()
        logs[ver].close()
        if proc.returncode != 0:
            failed.append(ver)
            print(f"[{ver}] FAILED (exit {proc.returncode}) — see run_{ver}.log")
        else:
            print(f"[{ver}] completed OK")
    if failed:
        sys.exit(f"Aborting: runs failed for versions {failed}")


def get_header(file):
    with open(file) as fh:
        header = fh.readline()
    return header.split("#")[1].split("\n")[0].split(" ")


def plot_RolligCombined(outfile='RolligCombined.pdf'):
    colors_mpl = plt.rcParams['axes.prop_cycle'].by_key()['color']
    RolligDir = "../tests/popsicle_semenov_photo_cr_full_RolligPDRBenchmark/RolligPDROutputs/"

    versions = ["V1", "V2", "V3", "V4"]
    col_titles = [
        r'$n = 10^3 \, \rm cm^{-3},\, \chi = 10 \; (V1)$',
        r'$n = 10^3 \, \rm cm^{-3},\, \chi = 10^5 \; (V2)$',
        r'$n = 10^{5.5} \, \rm cm^{-3},\, \chi = 10 \; (V3)$',
        r'$n = 10^{5.5} \, \rm cm^{-3},\, \chi = 10^5 \; (V4)$',
    ]
    panels = [
        (["Tgas", "Tdust"], [2, 3], r'$T_{\rm gas},\, T_{\rm dust} \, (\rm K)$',      True),
        (["H", "H2"],       [2, 3], r'$N_{\rm H},\, N_{\rm H_2} \, (\rm cm^{-2})$',   False),
        (["C+"],            [4],    r'$N_{\rm C^+} \, (\rm cm^{-2})$',                 False),
        (["C"],             [5],    r'$N_{\rm C} \, (\rm cm^{-2})$',                   False),
        (["CO"],            [6],    r'$N_{\rm CO} \, (\rm cm^{-2})$',                  False),
        (["E"],             [11],   r'$N_{\rm e^-} \, (\rm cm^{-2})$',                 False),
        (["H3+"],           [14],   r'$N_{\rm H_3^+} \, (\rm cm^{-2})$',               False),
        (["OH"],            [10],   r'$N_{\rm OH} \, (\rm cm^{-2})$',                  False),
        (["CH"],            [9],    r'$N_{\rm CH} \, (\rm cm^{-2})$',                  False),
        (["O"],             [7],    r'$N_{\rm O} \, (\rm cm^{-2})$',                   False),
    ]
    NROWS = len(panels)
    NCOLS = len(versions)

    fig, axs = plt.subplots(nrows=NROWS, ncols=NCOLS, figsize=(21.6, 43.2),
                            tight_layout=True, sharex=True)

    X_MIN, X_MAX = 1e-3, 10.0
    ymin = np.full((NROWS, NCOLS),  np.inf)
    ymax = np.full((NROWS, NCOLS), -np.inf)

    def _update_range(row, col, x, y):
        mask = (x >= X_MIN) & (x <= X_MAX) & np.isfinite(y) & (y > 0)
        if mask.any():
            ymin[row, col] = min(ymin[row, col], y[mask].min())
            ymax[row, col] = max(ymax[row, col], y[mask].max())

    for ci, version in enumerate(versions):
        krome_file = f"PDR_{version}"
        f = np.loadtxt(krome_file).T
        ntot = f[0]
        headers = get_header(krome_file)
        Av = f[4] * 6.289e-22          # column density → visual extinction
        dx = np.gradient(f[4]) / ntot  # path-length element dl = d(N_H) / n_tot

        for ri, (sp_list, _, _, raw) in enumerate(panels):
            for k, sp in enumerate(sp_list):
                color = colors_mpl[k]
                try:
                    ind = headers.index(sp)
                except ValueError:
                    continue
                # raw=True: temperature, plotted directly
                # raw=False: cumulative column density = cumsum(x_i * n_tot * dl)
                y = f[ind] if raw else np.cumsum(f[ind] * ntot * dx)
                axs[ri, ci].plot(Av, y, ls='-', lw=2.5, alpha=0.6,
                                 marker='o', markersize=8,
                                 markerfacecolor='none', markeredgecolor='k', mew=1.5,
                                 color=color)
                _update_range(ri, ci, Av, y)

        for rfile in glob.glob(RolligDir + f"*_{version}_Tgas.dat"):
            try:
                data = np.loadtxt(rfile, skiprows=1).T
            except Exception:
                continue
            for ri, (_, rcols, _, raw) in enumerate(panels):
                if not raw:
                    continue
                for k, cidx in enumerate(rcols):
                    if cidx >= data.shape[0]:
                        continue
                    axs[ri, ci].plot(data[1], data[cidx],
                                     color=colors_mpl[k], alpha=0.6, lw=1.0)

        for rfile in glob.glob(RolligDir + f"*_{version}.dat"):
            try:
                data = np.loadtxt(rfile, skiprows=1).T
            except Exception:
                continue
            dz = np.gradient(data[0])
            for ri, (_, rcols, _, raw) in enumerate(panels):
                if raw:
                    continue
                for k, cidx in enumerate(rcols):
                    if cidx >= data.shape[0]:
                        continue
                    axs[ri, ci].plot(data[1], np.cumsum(data[cidx] * dz),
                                     color=colors_mpl[k], alpha=0.6, lw=1.0)

        for ri in range(NROWS):
            axs[ri, ci].set_xscale('log')
            axs[ri, ci].set_yscale('log')
            axs[ri, ci].set_xlim(X_MIN, X_MAX)

    for ri in range(NROWS):
        for ci in range(NCOLS):
            if np.isfinite(ymin[ri, ci]) and ymax[ri, ci] > ymin[ri, ci]:
                axs[ri, ci].set_ylim(ymin[ri, ci], ymax[ri, ci])

    for ci, title in enumerate(col_titles):
        axs[0, ci].set_title(title, fontsize=24)

    for ri, (_, _, ylabel, _) in enumerate(panels):
        axs[ri, 0].set_ylabel(ylabel)

    for ci in range(NCOLS):
        axs[NROWS - 1, ci].set_xlabel(r'$A_{\rm V}$')

    fig.savefig(outfile, bbox_inches='tight')
    print(f"Saved {outfile}")


if __name__ == "__main__":
    os.chdir(BUILD_DIR)

    # Step 1: compile 4 versions sequentially
    for ver, chi0_val, ntot_val in VERSIONS:
        compile_version(ver, chi0_val, ntot_val)

    # Step 2: run all 4 in parallel
    procs, logs = run_all_parallel()

    # Step 3: wait and check
    wait_and_check(procs, logs)

    # Step 4: generate combined plot
    print("\nGenerating combined plot...")
    plot_RolligCombined()
