"""
Script to produce a plot of Tgas as a function of number density based on output from POPSICLE test
Author: Piyush Sharda (Leiden) 2024: sharda@strw.leidenuniv.nl
"""

import numpy as np
import matplotlib.pyplot as plt
import cmasher as cm

## Script to split the output files of the code into separate files based on metallicity
def split_by_metallicity(input_file, output_prefix,Zvals = [-2,-1,0]):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    print(Zvals)

    header = lines[0]  # <-- This is the added line to capture the header
    blocks = []
    current_block = []

    for line in lines:
        if line.strip() == "":
            if current_block:
                blocks.append(current_block)
                current_block = []
        else:
            current_block.append(line)
    
    # Add the final block if file doesn't end with a blank line
    if current_block:
        blocks.append(current_block)

    for i, block in enumerate(blocks):
        output_file = "{}_Z{}.dat".format(output_prefix, Zvals[i])
        with open(output_file, 'w') as f_out:
            f_out.writelines([header] + block)  # <-- This line ensures header is included
        print(f"Written {output_file} with {len(block)} lines")


Zvals = [-100,-6,-5,-4,-3,-2,-1,0]
Zlabels = [r'$-\infty$',r'$10^{-6} Z_{\odot}$',
           r'$10^{-5} Z_{\odot}$',r'$10^{-4} Z_{\odot}$',
           r'$10^{-3} Z_{\odot}$',r'$10^{-2} Z_{\odot}$',
           r'$10^{-1} Z_{\odot}$',r'$Z_{\odot}$']

# Example usage
split_by_metallicity("fort.22", "CollapseAB", Zvals=Zvals)
split_by_metallicity("fort.31", "CollapseCools", Zvals=Zvals)
split_by_metallicity("fort.911", "CollapseHeats", Zvals=Zvals)

#Now Plot
fig, axs = plt.subplots(ncols=1)
output_prefix = "CollapseAB"

for i, Zval in enumerate(Zvals):
    filename = "{}_Z{}.dat".format(output_prefix, Zval)
    f = np.loadtxt(filename).T
    axs.plot(f[0], f[2], label=Zlabels[i], color=cm.lavender(i/len(Zvals)), linewidth=3.0,alpha=0.6)

axs.set_xlabel(r'$n_{\rm H} \, (\rm cm^{-3})$', fontsize=20)
axs.set_ylabel(r'$T_{\rm gas} \ (\rm K)$', fontsize=20)
axs.set_yscale('log')
axs.set_xscale('log')

axs.legend(ncol=3)
fig.savefig("CollapseTgas_Metallicities", bbox_inches='tight')

