"""
Script to produce a plot of Tgas as a function of number density based on output from POPSICLE test
Author: Piyush Sharda (Leiden) 2024: sharda@strw.leidenuniv.nl
"""

import numpy as np
import matplotlib.pyplot as plt

defcolcycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

def read_file_with_breaks(filename, skip_header=0):
    with open(filename, 'r') as file:
        lines = file.readlines()

    data_sections = []
    current_section = []

    lines = lines[skip_header:]
    
    for line in lines:
        stripped_line = line.strip()
        if stripped_line:  # Non-empty line
            current_section.append(list(map(float, stripped_line.split())))
        else:  # Empty line
            if current_section:  # If there is data collected, save it to data_sections
                data_sections.append(np.array(current_section))
                current_section = []

    # Add the last section if not empty
    if current_section:
        data_sections.append(np.array(current_section))

    return data_sections


#Read in the output file
#Plot Tgas as a function of gas number density

filename = 'fort.22'
#Read in sections for different metallicity values
data_sections = read_file_with_breaks(filename, skip_header=1)

for i in range(0, len(data_sections)):
    plt.plot(np.log10(data_sections[i].T[0]), data_sections[i].T[1], c=defcolcycle[i], ls='solid')

plt.yscale('log')
plt.minorticks_on()
plt.ylim(1,4e3)
plt.ylabel('Tgas (K)')
plt.xlabel('log n (cm^-3)')
plt.grid()
plt.savefig('popsicle.png', bbox_inches='tight')

