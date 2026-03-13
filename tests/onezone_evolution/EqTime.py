import numpy as np
import matplotlib.pyplot as plt
import PhysicalConstantsCGS as const
import cmasher as cm
import re
from subprocess import run
from tabulate import tabulate  # for pretty output
import csv

colors_mpl = plt.rcParams["axes.prop_cycle"].by_key()["color"]

def get_Hnuclei(file):
    ab = open(file,'r')
    f = np.loadtxt(file).T
    header = ab.readline()
    headers = header.split("#")[1].split("\n")[0].split(" ")
    nHnuclei = 0.0
    cols = int(np.shape(f)[0])
    for i in range(2,cols):
        if(headers[i] == "H"):
            nHnuclei += f[i]
        elif(headers[i] == "H2"):
            nHnuclei += f[i] * 2.0
        elif(headers[i] == "H-"):
            nHnuclei += f[i]
        elif(headers[i] == "HD"):
            nHnuclei += f[i]
        elif(headers[i] == "OH"):
            nHnuclei += f[i]
        elif(headers[i] == "CH"):
            nHnuclei += f[i]
        elif(headers[i] == "CH2"):
            nHnuclei += f[i] * 2.0
        elif(headers[i] == "CH3"):
            nHnuclei += f[i] * 3.0
        elif(headers[i] == "CH4"):
            nHnuclei += f[i] * 4.0
        elif(headers[i] == "H20_total"):
            nHnuclei += f[i] * 2.0
        elif(headers[i] == "H+"):
            nHnuclei += f[i]
        elif(headers[i] == "H2+"):
            nHnuclei += f[i] * 2.0
        elif(headers[i] == "HD+"):
            nHnuclei += f[i]
        elif(headers[i] == "OH+"):
            nHnuclei += f[i]
        elif(headers[i] == "H20+"):
            nHnuclei += f[i] * 2.0
        elif(headers[i] == "H30+"):
            nHnuclei += f[i] * 3.0
        elif(headers[i] == "CH+"):
            nHnuclei += f[i]
        elif(headers[i] == "CH2+"):
            nHnuclei += f[i] * 2.0
        elif(headers[i] == "H3+"):
            nHnuclei += f[i] * 3.0

    ntot = f[0] #Get total number density since tabulated value scaled by this

    return nHnuclei*ntot

def get_Catoms(file):
    ab = open(file,'r')
    f = np.loadtxt(file).T
    header = ab.readline()
    headers = header.split("#")[1].split("\n")[0].split(" ")
    nCtot = 0.0
    cols = int(np.shape(f)[0])
    for i in range(2,cols):
        if(headers[i] == "C"):
            nCtot += f[i]
        elif(headers[i] == "CH"):
            nCtot += f[i]
        elif(headers[i] == "CH2"):
            nCtot += f[i]
        elif(headers[i] == "CH3"):
            nCtot += f[i]
        elif(headers[i] == "CH4"):
            nCtot += f[i]
        elif(headers[i] == "CO"):
            nCtot += f[i]
        elif(headers[i] == "CO2"):
            nCtot += f[i]
        elif(headers[i] == "C+"):
            nCtot += f[i]
        elif(headers[i] == "CO+"):
            nCtot += f[i]
        elif(headers[i] == "CH+"):
            nCtot += f[i]
        elif(headers[i] == "CH2+"):
            nCtot += f[i]

    ntot = f[0] #Get total number density since tabulated value scaled by this

    return nCtot*ntot

def get_Oatoms(file):
    ab = open(file,'r')
    f = np.loadtxt(file).T
    header = ab.readline()
    headers = header.split("#")[1].split("\n")[0].split(" ")
    nOtot = 0.0
    cols = int(np.shape(f)[0])
    for i in range(2,cols):
        if(headers[i] == "O"):
            nOtot += f[i]
        elif(headers[i] == "OH"):
            nOtot += f[i]
        elif(headers[i] == "H2O"):
            nOtot += f[i]
        elif(headers[i] == "O2"):
            nOtot += f[i] * 2.0
        elif(headers[i] == "CO"):
            nOtot += f[i]
        elif(headers[i] == "CO2"):
            nOtot += f[i] * 2.0
        elif(headers[i] == "O+"):
            nOtot += f[i]
        elif(headers[i] == "OH+"):
            nOtot += f[i]
        elif(headers[i] == "H2O+"):
            nOtot += f[i]
        elif(headers[i] == "H30+"):
            nOtot += f[i]
        elif(headers[i] == "O2+"):
            nOtot += f[i] * 2.0
        elif(headers[i] == "CO+"):
            nOtot += f[i]

    ntot = f[0] #Get total number density since tabulated value scaled by this

    return nOtot*ntot

def get_header(file):
    ab = open(file,'r')
    header = ab.readline()
    headers = header.split("#")[1].split("\n")[0].split(" ")
    return headers
        
           


def get_steady_state_time(t, x, tol=0.01, window=10):
    """
    Estimate the time when a quantity x(t) reaches steady-state.
    
    Parameters:
        t (ndarray): Time array (assumed sorted)
        x (ndarray): Abundance array (same shape as t)
        eps (float): Threshold for derivative (dx/dt) to be considered "steady"
        window (int): Number of consecutive steps to check for stability

    Returns:
        steady_time (float): Time when x(t) reaches steady-state
    """
    x_final = x[-1]
    threshold = tol * max(abs(x_final), 1e-30)

    for i in range(len(x)):
        if window == 0:
            if np.all(np.abs(x[i:] - x_final) < threshold):
                return t[i]
        else:
            if i + window <= len(x):
                if np.all(np.abs(x[i:i+window] - x_final) < threshold):
                    return t[i]
            #If we reach the end of the array without finding a steady state, return the last time
            else:
                return t[-1]
            
def get_teq_onezone(filename,verbose=True):
    f = np.loadtxt(filename).T
    headers = get_header(filename)
    nHnuclei = get_Hnuclei(filename)
    nCtot = get_Catoms(filename)
    nOtot = get_Oatoms(filename)

    filters = ["Tgas","H2","C","C+","CO"]
    column_indices = np.zeros(len(filters),dtype=int)
    #Read header and obtain the column number corresponding to the species
    for j,filter in enumerate(filters):
        ind = np.where(np.array(filter) == headers)[0][0]
        column_indices[j] = int(ind)

    #Now obtain the equilibrium time for the abundance
    tol = 0.1
    Tgas, nH2, nC, nCplus, nCO = f[column_indices[0]], f[column_indices[1]], f[column_indices[2]], f[column_indices[3]], f[column_indices[4]]
    teq_Temp, teq_H2, teq_C, teq_Cplus, teq_CO = get_steady_state_time(f[0], Tgas,tol=tol,window=100), get_steady_state_time(f[0], nH2,tol=tol,window=100), get_steady_state_time(f[0], nC,tol=tol,window=100), \
                                                get_steady_state_time(f[0], nCplus,tol=tol,window=100), get_steady_state_time(f[0], nCO,tol=tol,window=100)
    if(verbose is True):
        print("Equilibrium time for Tgas: {:.2f} Myr".format(teq_Temp/(1.e6*const.Year)))
        print("Equilibrium time for H2: {:.2f} Myr".format(teq_H2/(1.e6*const.Year)))
        print("Equilibrium time for C: {:.2f} Myr".format(teq_C/(1.e6*const.Year)))
        print("Equilibrium time for C+: {:.2f} Myr".format(teq_Cplus/(1.e6*const.Year)))
        print("Equilibrium time for CO: {:.2f} Myr".format(teq_CO/(1.e6*const.Year)))
    return teq_Temp, teq_H2, teq_C, teq_Cplus, teq_CO

def get_steadystate_values(filename):
    #Assume last time value is steady state value
    f = np.loadtxt(filename).T
    headers = get_header(filename)
    filters = ["Tgas","H2","C","C+","CO"]
    column_indices = np.zeros(len(filters),dtype=int)
    #Read header and obtain the column number corresponding to the species
    for j,filter in enumerate(filters):
        ind = np.where(np.array(filter) == headers)[0][0]
        column_indices[j] = int(ind)

    return f[column_indices[0]][-1], f[column_indices[1]][-1], f[column_indices[2]][-1], f[column_indices[3]][-1], f[column_indices[4]][-1]



def modify_krome_input_file(filepath, ntot_value, chi0_value, crate0_value):
    """
    Modify the Fortran file to update ntot, chi0, and crate_0 values.
    
    Args:
        filepath (str): Path to the Fortran file.
        ntot_value (float): Target gas number density.
        chi0_value (float): FUV scaling factor.
        crate0_value (float): Cosmic-ray ionization rate.
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Patterns and replacement logic
    def replace_line(pattern, replacement, lines):
        regex = re.compile(pattern)
        for i, line in enumerate(lines):
            if regex.search(line):
                lines[i] = replacement + "\n"
                break
        return lines

    lines = replace_line(r'^\s*ntot\s*=\s*[\deE.+-]+', f'    ntot = {ntot_value:.6e}', lines)
    lines = replace_line(r'^\s*chi0\s*=\s*[\deE.+-]+', f'  chi0 = {chi0_value:.6e}', lines)
    lines = replace_line(r'^\s*crate_0\s*=\s*[\deE.+-]+', f'  crate_0 = {crate0_value:.6e}', lines)

    with open(filepath, 'w') as f:
        f.writelines(lines)

    print(f"Updated {filepath} with ntot={ntot_value}, chi0={chi0_value}, crate_0={crate0_value}")

if __name__ == "__main__":

    # ntot_vals = np.logspace(0,7,num=100)
    # chi0_vals = np.logspace(0,5,num=10)  # FUV scaling factor
    # #crate_vals = np.logspace(-17,-13,num=1)  # Cosmic-ray ionization rate

    ntot_vals = np.array([1.e2,1.e3,1.e4,1.e5,1.e6,1.e7]) * 2.0
    #ntot_vals = np.array([1.e5]) * 2.0
    chi0_vals = [0,1,1e2,1e4]

    # Results list
    results = []

    #Loop over parameter space
    for ntot in ntot_vals:
        for chi0 in chi0_vals:
            # 1. Modify the input parameters
            crate = max(chi0*2e-16,2e-16)  # Ensure crate is at least 2e-16
            modify_krome_input_file("test.f90", ntot, chi0, crate) #Modify the input parameters
            #2. Compile the code
            compile_status = run(["make", "gfortran"], check=True)
            if compile_status.returncode != 0:
                print(f"Compilation failed for ntot={ntot}, chi0={chi0}, crate={crate}")
                continue
            # 3. Run the executable
            run_status = run(["./test"])
            if run_status.returncode != 0:
                print(f"Run failed for ntot={ntot}, chi0={chi0}, crate={crate}")
                continue
            # 4. Process the output files
            #split_by_metallicity("fort.22", "OneZoneAB")
            filename = "fort.22"
            #5. Extract the equilibrium time for the species of interest
            teq_Temp, teq_H2, teq_C, teq_Cplus, teq_CO = get_teq_onezone(filename)
            #6. Store the equilibrium values
            eq_value_Temp, eq_value_H2, eq_value_C, eq_value_Cplus, eq_value_CO = get_steadystate_values(filename)
            # 5. Store result
            results.append({
                "ntot": ntot,
                "chi0": chi0,
                "crate": crate,
                "Teq_Temp": teq_Temp,
                "Teq_H2": teq_H2,
                "Teq_C": teq_C,
                "Teq_C+": teq_Cplus,
                "Teq_CO": teq_CO,
                "Eq_Temp": eq_value_Temp,
                "Eq_H2": eq_value_H2,
                "Eq_C": eq_value_C,
                "Eq_Cplus": eq_value_Cplus,
                "Eq_CO": eq_value_CO
            })



    # Print the results as a nice table
    headers = ["ntot", "chi0", "crate", "Teq_Temp", "Teq_H2", "Teq_C", "Teq_C+", "Teq_CO"]
    table = [[
        r["ntot"], r["chi0"], r["crate"],
        f"{r['Teq_Temp']:.3e}", f"{r['Teq_H2']:.3e}", f"{r['Teq_C']:.3e}",
        f"{r['Teq_C+']:.3e}", f"{r['Teq_CO']:.3e}"
    ] for r in results]

    print(tabulate(table, headers=headers, tablefmt="grid"))

    # write to CSV
    with open("equilibrium_times.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)