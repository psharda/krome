import numpy as np
import matplotlib.pyplot as plt
import glob
import os

# ── helpers ──────────────────────────────────────────────────────────────────

def get_header(file):
    with open(file) as fh:
        header = fh.readline()
    return header.split("#")[1].split("\n")[0].split(" ")


# ── load KROME output ────────────────────────────────────────────────────────

krome_file = os.path.join(os.path.dirname(__file__), "PDR_Z0")
f = np.loadtxt(krome_file).T
ntot = f[0]
headers = get_header(krome_file)

N_to_Av = 6.289e-22
Av = f[4] * N_to_Av


def col_index(name):
    return np.where(np.array(name) == headers)[0][0]


# ── 4×2 figure ───────────────────────────────────────────────────────────────

fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(21.6, 9.6),
                        tight_layout=True, sharex=True)

# Panel layout (row, col) → (krome species, Rollig col index, y-label)
panels = [
    # top row
    (0, 0, ["H", "H2"], [3, 2],   r'$n_{\rm H},\,n_{\rm H_2}\,(\rm cm^{-3})$'),
    (0, 1, ["C+"],       [4],     r'$n_{\rm C^+}\,(\rm cm^{-3})$'),
    (0, 2, ["C"],        [5],     r'$n_{\rm C}\,(\rm cm^{-3})$'),
    (0, 3, ["CO"],       [6],     r'$n_{\rm CO}\,(\rm cm^{-3})$'),
    # bottom row
    (1, 0, ["H3+"],      [14],    r'$n_{\rm H_3^+}\,(\rm cm^{-3})$'),
    (1, 1, ["E"],        [11],    r'$n_{\rm e}\,(\rm cm^{-3})$'),
    (1, 2, ["CH"],       [9],     r'$n_{\rm CH}\,(\rm cm^{-3})$'),
    (1, 3, ["OH"],       [10],    r'$n_{\rm OH}\,(\rm cm^{-3})$'),
]

# Plot KROME data
for row, col, species_list, _, ylabel in panels:
    ax = axs[row, col]
    for species in species_list:
        idx = col_index(species)
        ax.plot(Av, f[idx] * ntot, ls='none', marker='o', markersize=4,
                alpha=0.7, label=species)
    ax.set_ylabel(ylabel)

# Rollig benchmark (black lines)
#rollig_dir = "/Users/smenon/Desktop/Work/PhotoMetalImplement/Tests/Full_PDR/RolligPDR/"
rollig_dir = "../tests/popsicle_semenov_photo_cr_full_RolligPDRBenchmark/RolligPDROutputs/"
for rollig_file in sorted(glob.glob(rollig_dir + "*_V1.dat")):
    try:
        d = np.loadtxt(rollig_file, skiprows=1).T
    except Exception as e:
        print(f"skip {os.path.basename(rollig_file)}: {e}")
        continue
    kw = dict(c='k', alpha=0.5, lw=1.0)
    for row, col, _, rollig_cols, _ in panels:
        ax = axs[row, col]
        try:
            for rc in rollig_cols:
                ax.plot(d[1], d[rc], **kw)
        except IndexError:
            pass

# Axes formatting
for ax in axs.flat:
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1.e-3, 10.0)

# y-limits for first four panels
axs[0, 0].set_ylim(1.e-2, 2.e6)
axs[0, 1].set_ylim(1.e-6, 1.e3)
axs[0, 2].set_ylim(1.e-6, 1.e3)
axs[0, 3].set_ylim(1.e-6, 1.e3)

# x-labels on bottom row
for col in range(4):
    axs[1, col].set_xlabel(r'$A_{\rm V}$')

# legend on H/H2 panel
axs[0, 0].legend(fontsize=9)

fig.savefig(os.path.join(os.path.dirname(__file__), "RolligPDR.pdf"),
            bbox_inches='tight')
print("Saved RolligPDR.pdf")
