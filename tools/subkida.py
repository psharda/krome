#
# THIS SCRIPT CONVERT A KIDA DATABASE INTO A KROME-READABLE FILE.
# THIS FILE IS PROVIDED AS IT IS, PLEASE CAREFULLY CHECK THE OUTPUT PRODUCED.
# SEE THE DISCLAIMER AT THE END OF THIS SECTION.
#
# This script is a part of KROME.
#
# see https://bitbucket.org/tgrassi/krome
#
# Contact Tommaso Grassi for further details
# tgrassi@usm.lmu.de


# convert from KIDA file to KROME network (including a set of option to chose a subset)

# input filename (in KIDA format, see http://kida.obs.u-bordeaux1.fr/networks.html)
fname = "kida_demo.dat"
# output filename for KROME
foutname = "react_subkida"
# filename for multiple reactions (same reactants and same products)
fmultname = "react_multi"

# several options to select your own reaction subset
avoid = []  # avoid atoms (can be empty)
# use atoms (can be empty, e.g. ["H", "He"]. "e"=electron, "grain", "_ortho", "_para", ...)
# example use = ["H", "He", "D", "e", "_para", "_ortho", "_meta", "GRAIN"]
use = []
maxatoms = 999  # maximum number of atoms
ions = True  # use ions
anions = True  # use anions
Tmin = -1e99  # minimum temperature
Tmax = 1e99  # maximum temperature
include = []  # include these species (ignored if empty)
exclude = []  # species to exclude (ignored if empty)
excludein = []  # exclude species that contain a specific case-sensitive string (ignored if empty e.g. l- for linear mols)
multiple = True  # include multiple reactions (same reactants and same products)
extend_single = True  # remove temperature limits for reactions with Tmin=Tmax

# Recommendation is the recommendation given by experts in KIDA. (from KIDA help)
# 0 means that the value is not recommended.
# 1 means that there is no recommendation (experts have not looked at the data).
# 2 means that the value has been validated by experts over the Trange
# 3 means that it is the recommended value over Trange.
recom = [1, 2, 3]  # recomandations to include (see above)

# Processes (selected using fitting formula)
# 0: dust reaction as in Majumdar+2017
# 1: Cosmic-ray ionization (direct and undirect)
# 2: Photo-dissociation (Draine)
# 3: Kooij
# 4: ionpol1
# 5: ionpol2
processes = [0, 1, 2, 3, 4, 5]  # processes included (see above)

# variables for cosmic rays and photochemistry
CRvar = "user_crflux"  # name of the CR flux variable
Avvar = "user_Av"  # name of the Av variable

# ******************************************************************
# ******************************************************************
# ******************************************************************
#  MODIFY BELOW THIS LINE ONLY IF YOU KNOW WHAT TO DO!
# ******************************************************************
# ******************************************************************
# ******************************************************************

use = [x.upper() for x in use]
avoid = [x.upper() for x in avoid]
exclude = [x.upper() for x in exclude]
include = [x.upper() for x in include]
recom = [int(x) for x in recom]
processes = [int(x) for x in processes]

import sys
import datetime

# define number of reactants and products
reactants_number = 3
products_number = 5


# ************************
# function to check if argument is a number
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


# *********************
# format items list as tlen-sized columns
def tabrow(items, tlen=20):
    frow = ""
    for item in items:
        frow += str(item) + " "*(tlen-len(str(item)))
    return frow


# ************************
# molecule class
class Mol:
    name = ""  # name
    atoms = []  # list of atoms
    natoms = 0  # count atoms
    charge = 0  # charge
    emol = []  # exploded molecule HCO2H = [H,H,C,O,O]

    # method to parse the name
    def parse(self):
        rms = ["c-", "l-", "(1)", "+", "-", "(2D)", "(1D)"]
        adic = [str(jj) for jj in range(30)]  # numbers
        adic += ["E", "H", "D", "HE", "LI", "C", "N", "O", "F", "NE", "NA", "MG", "AL",
                 "SI", "P", "S", "CL", "CA", "K", "MN", "FE", "NI"]  # elements
        adic += ["_META", "_ORTHO", "_PARA", "GRAIN"]
        adic = sorted(adic, key=lambda xa: len(xa), reverse=True)  # sort by length
        for kop, vop in {"m": "_meta", "o": "_ortho", "p": "_para"}.iteritems():
            if self.name.startswith(kop):
                self.name = self.name[1:] + vop
        name = self.name  # copy name locally
        cname = name
        for xx in rms:
            cname = cname.replace(xx, "")
        cname = cname.upper()  # copy name and remove charge
        amol = []
        atoms = []
        # loop to find multiple objects
        for j in range(10):
            # loop on atoms
            for a in adic:
                # serach for element
                if a in cname:
                    amol.append([cname.index(a), a])  # append position and name
                    cname = cname.replace(a, "X"*len(a), 1)  # replace element with XX
                    if (a not in atoms) and (not is_number(a)):
                        atoms.append(a)  # store atom
        amol = sorted(amol, key=lambda xa: xa[0])  # sort by position
        emol = []  # exploded molecule
        if cname.replace("X", "").strip() != "":
            print "ERROR: problem when parsing " + self.name
            print name, cname
            return None

        ii = 0
        while True:
            mult = 1  # multeplicity (e.g. H2 => mult=2)
            if ii >= len(amol):
                break  # end loop
            # if number skip
            if is_number(amol[ii][1]):
                ii += 1
                continue
            # if next exists
            if ii < len(amol)-1:
                # if next is number => is multeplicity
                if is_number(amol[ii+1][1]):
                    mult = int(amol[ii+1][1])
            emol += ([amol[ii][1]] * mult)  # add to exploded
            ii += 1
        self.emol = emol  # exploded
        self.atoms = atoms  # atoms
        self.natoms = len(emol)  # count atoms
        self.charge = name.count("+") - name.count("-")  # charge

        # returns self.name if correctly parsed, otherwise returns None
        return self.name

# **************************************
# **************************************
# **************************************
# Format : 3(a11) 1x 5(a11) 1x 3(e10.3 1x) 2(e8.2) 1x a4 i3 2(i7) i3 4x 2(i2) 1x i2
# Reactants   Products  alpha  beta  gamma  F g Type_of_uncertainty   itype    Tmin  Tmax  Formula   Number     Number_of_(alpha, beta, gamma)
# Recommendation


print "******************************"
print "****** RUNNING SUB KIDA ******"
print "******************************"

print "reading file "+fname+", wait..."

okcount = totcount = trangecount = single_count = 0
nhist = dict()
idxdic = dict()
trange = dict()
tsingle = dict()
formulahist = dict()
fh = open(fname, "rb")
fout = open(foutname, "w")
fmult = open(fmultname, "w")
fmt = [11]*reactants_number + [1] + products_number*[11] + [1] + 3*[11] \
      + [8, 9] + [1, 4, 3] + 2*[7] + [3, 6, 3, 2]
keys = ["R"+str(i) for i in range(reactants_number)] + ["x"] + \
        ["P"+str(i) for i in range(products_number)] + ["x"] + ["a", "b", "c"] \
        + ["F", "g"] + ["x", "unc", "type"] + ["tmin", "tmax"] \
        + ["formula", "num", "subnum", "recom"]
reacts = []
fout.write("# reaction subset from KIDA created with KROME (see Grassi+2014)\n")

# time as string
date_fmt = "%Y-%m-%d %H:%M:%S"
date_string = "# creation date:" + datetime.datetime.now().strftime(date_fmt) + "\n"

# write creation time to output files
fout.write(date_string)
fmult.write(date_string)

if len(avoid) > 0:
    fout.write("#avoid=" + (",".join(avoid)) + "\n")
if len(use) > 0:
    fout.write("#use=" + (",".join(use)) + "\n")
if len(exclude) > 0:
    fout.write("#exclude=" + (",".join(exclude)) + "\n")

fout.write("#recom=" + (",".join([str(x) for x in recom])) + "\n")
fout.write("#processes=" + (",".join([str(x) for x in processes])) + "\n")

if maxatoms < 100:
    fout.write("#maxatoms=" + str(maxatoms) + "\n")
if not ions:
    fout.write("#no ions\n")
if not anions:
    fout.write("#no anions\n")
fout.write("#Tmin=" + str(Tmin) + " Tmax=" + str(Tmax) + "\n")
fout.write("@common:" + CRvar + "," + Avvar + "\n")
fout.write("@format:idx,R,R,R,P,P,P,P,P,Tmin,Tmax,rate\n")

formula_not_found_count = 0
rems = ["Photon", "CRP", "CRPHOT", "CR", ""]
for row in fh:
    srow = row.strip()
    if srow == "":
        continue  # skip empty lines
    if srow.startswith("#"):
        continue  # skip comments
    if srow.startswith("!"):
        continue  # skip comments
    totcount += 1

    p = 0
    arow = dict()
    for i in range(len(fmt)):
        arow[keys[i]] = srow[p:p+fmt[i]].strip()
        p += fmt[i]
    RR = ",".join([arow["R"+str(i)] for i in range(reactants_number)])
    PP = ",".join([arow["P"+str(i)] for i in range(products_number)])
    TT = arow["tmin"] + "," + arow["tmax"]

    # Formula is a number that referes to the formula needed to compute the rate
    # coefficient of the reaction.
    # see http://kida.obs.u-bordeaux1.fr/help
    # 1: Cosmic-ray ionization (direct and undirect)
    # 2: Photo-dissociation (Draine)
    # 3: Kooij
    # 4: ionpol1
    # 5: ionpol2
    # 6: Troe fall-off (NOT SUPPORTED!)
    arow["formula"] = int(arow["formula"])
    KK = None
    if arow["formula"] == 1:
        KK = arow["a"]+"*"+CRvar
    elif arow["formula"] == 2:
        KK = arow["a"]
        if float(arow["c"]) != 0e0:
                KK += "*exp(-"+arow["c"]+"*"+Avvar+")"
    elif arow["formula"] in [0, 3]:
        if arow["formula"] == 0:
            print "WARNING: found rate with Formula type 0 treated as Kojii (i.e. type 3)"
        KK = arow["a"]
        if float(arow["b"]) != 0e0:
            KK += "*(T32)**("+arow["b"]+")"
        if float(arow["c"]) != 0e0:
            KK += "*exp(-"+arow["c"]+"*invT)"
    elif arow["formula"] == 4:
        KK = arow["a"]
        if float(arow["b"]) != 1e0:
            KK += "*"+arow["b"]
        gpart = ""
        if float(arow["c"]) != 0e0:
            gpart = "+ 0.4767d0*"+arow["c"]+"*sqrt(3d2*invT)"
        KK += "*(0.62d0 "+gpart+")"
    elif arow["formula"] == 5:
        KK = arow["a"]
        if float(arow["b"]) != 1e0:
            KK += "*"+arow["b"]
        gpart = ""
        if float(arow["c"]) != 0e0:
            gpart = "+ 0.0967d0*"+arow["c"]+"*sqrt(3d2*invT) + "
            gpart += arow["c"]+"**2*28.501d0*invT"
            KK += "*(1d0 "+gpart+")"
    else:
        print "WARNING: Formula not found!", arow["formula"]
        formula_not_found_count += 1
        print srow[:60] + "..."
        continue

    KK = KK.replace("--", "+").replace("++", "+").replace("-+", "-").replace("+-", "-")

    ok = True
    if int(arow["formula"]) not in processes:
        ok = False
        print "skip reaction with process", int(arow["formula"])
    if int(arow["recom"]) not in recom:
        print "skip reaction according to recom", int(arow["recom"])
        ok = False
    if (arow["formula"] != 1) and (arow["formula"] != 2):
        if float(arow["tmin"]) < Tmin:
            ok = False
            print "skip reaction according to Tmin", Tmin
        if float(arow["tmax"]) > Tmax:
            ok = False
            print "skip reaction according to Tmax", Tmax
    if not ok:
        continue

    RR_obj = []
    for i in range(reactants_number):
        if not ok:
            break
        mymol = Mol()
        R = arow["R"+str(i)]
        if R in rems:
            continue
        if R in exclude:
            ok = False
        # check if string from stringlist excludein is in R
        if len(excludein) > 0:
            for ex in excludein:
                if ex in R:
                    ok = False
                    break

        mymol.name = R

        check = mymol.parse()
        if check is None:
            print "ERROR: problem when parsing a species in line"
            print srow
            sys.exit()

	if include:
		if mymol.name not in include:
			ok = False

        RR_obj.append(mymol)

        if mymol.natoms > maxatoms:
            ok = False
        if ("+" in R) and not ions:
            ok = False
        if ("-" in R) and not anions:
            ok = False
        for x in avoid:
            if x in mymol.atoms:
                ok = False
        if len(use) > 0:
            for x in mymol.atoms:
                if x not in use:
                    ok = False

    PP_obj = []
    for i in range(products_number):
        if not ok:
            break
        mymol = Mol()
        P = arow["P"+str(i)]
        if P in rems:
            continue
        if P in exclude:
            ok = False
        # check if string from stringlist excludein is in P
        if len(excludein) > 0:
            for ex in excludein:
                if ex in P:
                    ok = False
                    break
        mymol.name = P
        check = mymol.parse()

        if check is None:
            print "ERROR: problem when parsing a species in line"
            print srow
            sys.exit()

	if include:
		if mymol.name not in include:
			ok = False

        PP_obj.append(mymol)

        if mymol.natoms > maxatoms:
            ok = False
        if ("+" in P) and not ions:
            ok = False
        if ("-" in P) and not anions:
            ok = False
        for x in avoid:
            if x in mymol.atoms:
                ok = False
        if len(use) > 0:
            for x in mymol.atoms:
                if x not in use:
                    ok = False
    if not ok:
        continue

    # check reactions with Tmax==Tmin
    if float(arow["tmin"]) == float(arow["tmax"]):
        single_count += 1
        if extend_single:
            arow["tmin"] = 0e0
            arow["tmax"] = 1e10
        else:
            if float(arow["tmin"]) not in tsingle:
                tsingle[float(arow["tmin"])] = 0
            tsingle[float(arow["tmin"])] += 1

    if arow["formula"] not in formulahist:
        formulahist[arow["formula"]] = 0
    formulahist[arow["formula"]] += 1

    isMult = False
    if arow["num"] in nhist:
        sameTRange = False
        if (float(arow["tmin"]) > trange[arow["num"]][0]) \
                and (float(arow["tmin"]) < trange[arow["num"]][1]):
            sameTRange = True
        if (float(arow["tmax"]) > trange[arow["num"]][0]) \
                and (float(arow["tmax"]) < trange[arow["num"]][1]):
            sameTRange = True
        if float(arow["tmin"]) == trange[arow["num"]][0]:
            sameTRange = True
        if float(arow["tmax"]) == trange[arow["num"]][1]:
            sameTRange = True
        if not sameTRange:
            trangecount += 1
        if (not multiple) and sameTRange:
            continue
        if sameTRange:
            nhist[arow["num"]].append(arow["subnum"])
            idxdic[arow["num"]] = okcount + 1
            isMult = True
    else:
        nhist[arow["num"]] = [arow["subnum"]]
        trange[arow["num"]] = [float(arow["tmin"]), float(arow["tmax"])]

    RR_names = ",".join([x.name for x in RR_obj])
    PP_names = ",".join([x.name for x in PP_obj])
    if len(RR_obj) < reactants_number:
        RR_names += ","*(reactants_number - len(RR_obj))
    if len(PP_obj) < products_number:
        PP_names += ","*(products_number - len(PP_obj))

    krow = (RR_names + "," + PP_names + "," + TT + "," + KK + "\n")
    for x in rems:
        krow = krow.replace(x, "")

    # @format:idx,R,R,R,P,P,P,P,P,Tmin,Tmax,rate
    akrow = krow.split(",")
    akrow_new = []
    for prt in akrow:
        if len(prt) > 2:
            if (prt[:2] == "c-") or (prt[:2] == "l-"):
                lprt = list(prt)
                lprt[1] = "_"
                prt = "".join(lprt)
        akrow_new.append(prt)
    krow = (",".join(akrow_new))

    okcount += 1
    if isMult:
        fmult.write(str(okcount) + "," + krow)
    fout.write(str(okcount) + "," + krow)
    reacts.append(arow)

fmult.close()
fout.close()
fh.close()


# FINAL OUTPUT
if okcount == 0:
    print "********** WARNING! **********"
    print "No reactions matching your criteria!!!"
    sys.exit()

multi = [str(k) + "("+str(idxdic[k])+") [" +
         (",".join(v)) + "]" for (k, v) in nhist.iteritems() if len(v) > 1]

if len(multi) > 1:
    print "WARNING: there are " + str(len(multi)) + " reactions with multiple values"
    print " Reactions with multiple values written in", fmultname

print "Total reactions:", totcount
print "Rections INCLUDED:", okcount
print "Rections NOT INCLUDED:", totcount - okcount
print "Multiple reactions (same reactants and products):", len(multi)
print "Formula not found reactions:", formula_not_found_count
print "Joined Trange reactions:", trangecount
if len(tsingle) > 0:
    print "WARNING: Found reactions with Tmin==Tmax:", sum(tsingle.values()), "as"
    print "----------------------------------------"
    print " " + tabrow(["Tmin=Tmax", "count"])
    print "----------------------------------------"
    for k, v in tsingle.iteritems():
        print " " + tabrow([str(k), v])
    print "----------------------------------------"
if extend_single and (single_count > 0):
    print "WARNING: temperature limits removed for " + str(single_count) \
          + " reactions with Tmin=Tmax"
print "Formula count per type:"

rtype = {0: "Dust/Special",
         1: "CR ioniz",
         2: "Photo-diss",
         3: "Kooij",
         4: "ionpol1",
         5: "ionpol2"}
if len(formulahist) > 0:
    for k, v in formulahist.iteritems():
        print " " + tabrow([rtype[k], ":", v])
print "File written in:", foutname
print "Reactions with multiple values written in", fmultname

print "Everything done! Bye!"
