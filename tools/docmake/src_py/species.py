import utils, itertools, os, sys


class species():
    name = nameLatex = nameFile = nameHtml = nameHref = None

    # *****************
    # constructor
    def __init__(self, speciesName, atomSet, skipPhotochemistry=False):

        # check for upper case atoms, e.g. HE instead of He
        # only if species (without signs) is not in the list
        if speciesName.replace("+", "").replace("-", "") not in atomSet.keys():
            # loop on "atoms" names
            for atom in atomSet.keys():
                # if found lowercase replace it
                if (atom.upper() != atom) and (atom.upper() in speciesName):
                    orgSpeciesName = speciesName
                    speciesName = speciesName.replace(atom.upper(), atom)
                    print("WARNING: " + atom + " found upper case in " + orgSpeciesName \
                          + "! Replaced as " + speciesName)

        if speciesName.startswith("J"):
            speciesName = speciesName[1:] + "_dust"

        # replace linear and chain prefix with underscore
        isoreplace = {"c-": "c_", "l-": "l_"}
        for (iso, isorep) in isoreplace.items():
            if speciesName.startswith(iso) and len(speciesName) > len(iso):
                speciesName = speciesName.replace(iso, isorep)

        self.name = speciesName

        self.nameLatex = self.getLatexName()
        self.nameFile = self.nameHash = self.name.lower().replace("+", "j").replace("-", "w")
        self.nameHtml = self.getHtmlName()
        self.nameHref = self.getHrefName()

        me = 9.10938356e-28  # electron mass, g

        # store keys sorted by inverse length
        atoms = sorted(atomSet.keys(), key=lambda x: len(x), reverse=True)
        # produce unique character combinations
        alpha = ["".join(x) for x in list(itertools.product("XYZ", repeat=4))]
        # check to have enough combinations
        if len(atoms) > len(alpha):
            sys.exit("ERROR: in species parser alpha needs to be extended!")

        # local species name
        specName = speciesName  # .upper()

        # compute charge
        positiveCount = specName.count("+")
        negativeCount = specName.count("-")
        self.charge = positiveCount - negativeCount
        if speciesName.upper() == "E":
            self.charge = -1
        # replace signs
        specName = specName.replace("+", "").replace("-", "")

        # replace atoms slash*separated
        for i in range(len(atoms)):
            specName = specName.replace(atoms[i], "/" + alpha[i] + "/")
        # replace double slashes
        while "//" in specName:
            specName = specName.replace("//", "/")
        # split at slashes
        aspec = [x for x in specName.split("/") if x != ""]

        # search for number and when found multiplicate previous non-number
        exploded = []
        for a in aspec:
            if utils.isNumber(a):
                for j in range(int(a) - 1):
                    exploded.append(aold)
            else:
                exploded.append(a)
            aold = a

        # store exploded with real atom names
        try:
            self.exploded = [atoms[alpha.index(x)] for x in exploded]
        except:
            print("ERROR: wanted to parse ", self.name)
            print(" but something went wrong with ", exploded)
            print(" Available atoms are:", atoms)
            print(" Add to atomlist.dat if needed.")
            sys.exit()
        self.explodedFull = self.exploded + (["+"] * positiveCount) + (["-"] * negativeCount)

        # store atoms
        nonAtoms = ["+", "-"]
        self.atoms = [x for x in self.exploded if (not (x in nonAtoms) and not (utils.isNumber(x)))]

        # compute mass using dictionary as reference
        self.mass = sum([atomSet[x] for x in self.exploded])
        if speciesName.upper() != "E": self.mass -= self.charge * me

        # store radius
        self.radius = self.getRadius(self.name)
        # skip photochemistry if required
        if not skipPhotochemistry:
            # init photochem variables
            self.xsecs = dict()
            self.phrates = dict()
            self.phrates_extinction = dict()

            # load xsec from file
            self.loadXsecLeiden()

            # compute photo rates
            self.computePhRates()

    # **********************
    # load xsecs, structure is a dict
    # database_name
    #   |---energy---[list of energies, eV]
    #   |---phd/phi---[list of xsecs, eV]
    def loadXsecLeiden(self):

        fname = "xsecs/" + self.name + ".dat"
        if not os.path.exists(fname): return

        # get file size in MB
        sizeMB = round(os.path.getsize(fname) / 1024 ** 2, 1)

        # check if xsec size is relatively large
        if sizeMB > 10.0:
            print("WARNING: " + fname + " is quite large (" + str(sizeMB) + " MB)")
            # ask if to load xsec, N is default
            while True:
                reply = input("Load Xsec file? (y/N) ").lower().strip()
                if reply == "y":
                    break
                if reply == "n" or reply == "":
                    return

        self.xsecs["leiden"] = dict()

        clight = 2.99792458e10  # cm/s
        hplanck = 4.135667662e-15  # eV*s

        # loop on file rows
        for row in open(fname):
            srow = row.strip()
            if srow == "": continue
            # look for format
            if srow.startswith("# wavelength"):
                fmt = [x for x in srow.split(" ") if x != ""][1:]
                # init arrays
                for fm in fmt:
                    self.xsecs["leiden"][fm] = []
            if srow.startswith("#"): continue
            arow = [float(x) for x in srow.split(" ") if (x != "")]
            # loop on data
            for i in range(len(fmt)):
                self.xsecs["leiden"][fmt[i]].append(arow[i])

        self.xsecs["leiden"]["energy"] = [clight * hplanck / (wl * 1e-7) for wl in self.xsecs["leiden"]["wavelength"]]

        # reverse data
        for k, v in self.xsecs["leiden"].items():
            self.xsecs["leiden"][k] = v[::-1]

    # **********************
    # compute photo rates, structure is a dict
    # database_name
    #   |---phd/phi
    #      |---radiation_fields---(rate,1/s)
    def computePhRates(self):
        import numpy as np
        from math import pi
        from photo import Jdraine, JBB, intJdraine, intJBB

        hplanck = 4.135667662e-15  # eV*s

        # radiation types
        frads = {"Draine1978": Jdraine,
                 "BB@4e3K": JBB,
                 "BB@1e4K": JBB,
                 "BB@2e4K": JBB}

        if len(self.xsecs) != 0:
            print("  Computing {0} photochemical rates and extinction coefficients...".format(self.name))

        # loop on databases
        for (db, data) in self.xsecs.items():
            self.phrates[db] = dict()
            xdata = data["energy"]
            # loop on rates
            for k, ydata in data.items():
                # skip these keys
                if k in ["energy", "wavelength", "photoabsorption"]: continue
                if len(ydata) == 0: continue
                self.phrates[db][k] = dict()
                # loop on radiation fluxes
                for radName, Jfun in frads.items():
                    if radName.startswith("BB@"):
                        Tbb = float(radName.replace("BB@", "").replace("K", ""))
                        fscale = intJdraine() / intJBB(Tbb)
                    fdata = []
                    # loop on energy range
                    for i in range(len(xdata)):
                        sigma = ydata[i]  # cm2
                        energy = xdata[i]  # eV
                        if radName.startswith("BB@"):
                            Jrad = fscale * Jfun(energy, Tbb)  # eV/cm2/sr
                        else:
                            Jrad = Jfun(energy)  # eV/cm2/sr
                        fdata.append(Jrad * sigma / energy)

                    # compute integral and rate, 1/s
                    kph = 4e0 * pi * np.trapz(fdata, xdata) / hplanck

                    # store rate, 1/s
                    self.phrates[db][k][radName] = kph

        # compute b in a*exp(-b*Av) WITHOUT SCATTERING!
        self.computeExtinction()

    # *******************
    # compute rate extinction WITHOUT SCATTERING!
    def computeExtinction(self):
        import numpy as np
        from scipy.interpolate import interp1d
        from scipy.optimize import curve_fit
        from math import pi, exp
        from photo import Jdraine, JBB, intJdraine, intJBB

        # file with dust information, from Draine
        # https://www.astro.princeton.edu/~draine/dust/dustmix.html
        fname = "kext_albedo_WD_MW_3.1_60_D03.all"

        # some constants
        clight = 2.99792458e10  # speed of light, cm/s
        hplanck = 4.135667662e-15  # planck constant, eV*s
        mp = 1.6726219e-24  # proton mass, g

        energyList = []
        kabsList = []
        cos2List = []
        albedoList = []
        # flag to skip header
        inRead = False
        # loop to read dust file
        for row in open(fname):
            srow = row.strip()
            # starts to read when this string has found
            if srow.startswith("-----------"):
                inRead = True
                continue
            # skip if still inside comments
            if not inRead:
                continue

            # read data
            arow = [x for x in srow.split(" ") if x != ""]
            # convert data to float, comments are not included
            (wavelength, albedo, cosAvg, cext, kabs, cos2Avg) = [float(x) for x in arow[:6]]

            # convert wavelenght into energy, note wavelenght is descending
            # hence no need to reverse array after conversion
            energyList.append(clight * hplanck / (wavelength * 1e-4))
            # store opacity, cm2/g
            kabsList.append(kabs)
            # store averaged scattering angle squared <cos2(theta)>
            cos2List.append(cos2Avg)
            # store albedo
            albedoList.append(albedo)

        # interpolate stored quantities as functions of energy
        fkabs = interp1d(energyList, kabsList)
        fcos2 = interp1d(energyList, cos2List)
        falbedo = interp1d(energyList, albedoList)

        # set a gas model to propagate radiation
        ngas = 1e0  # gas density, 1/cm3
        d2g = 0.00806451612  # dust/gas mass ratio, 1/124
        mu = 1e0  # mean mol weight
        rho_dust = ngas * mp * mu * d2g  # dust density, cm2/g
        minAv = 0.1  # min Av of the seminfinite slab
        maxAv = 3e0  # max Av of the seminfinite slab

        # Av->col, cm2
        av2ncol = 1.6e21

        # compute min/max position
        xmin = minAv * av2ncol / ngas
        xmax = maxAv * av2ncol / ngas

        # number of grid points in space
        imax = 30

        # loop on database
        for db, data in self.xsecs.items():
            self.phrates_extinction[db] = dict()
            xdata = data["energy"]
            # loop on rates
            for (k, ydata) in data.items():
                # skip these keys
                if k in ["energy", "wavelength", "photoabsorption"]:
                    continue
                # skip empty data
                if len(ydata) == 0:
                    continue

                self.phrates_extinction[db][k] = dict()

                kphList = []
                avList = []

                # loop on grid points
                for xpos in [i * (xmax - xmin) / (imax - 1) + xmin for i in range(imax)]:
                    # loop on radiation fluxes
                    edata = []
                    fdata = []
                    # loop on energy range
                    for i in range(len(xdata)):
                        sigma = ydata[i]  # cm2
                        energy = xdata[i]  # eV
                        if energy < min(energyList) or energy > max(energyList): continue
                        edata.append(energy)
                        tau = rho_dust * xpos * fkabs(energy)
                        # Jrad = fscale*JBB(energy, 4e3)*exp(-tau)*fcos2(energy)/(1e0-falbedo(energy))
                        Jrad = Jdraine(energy) * exp(-tau)
                        fdata.append(Jrad * sigma / energy)

                    # compute integral and rate, 1/s
                    kph = 4e0 * pi * np.trapz(fdata, edata) / hplanck
                    kphList.append(kph)
                    avList.append(xpos * ngas / av2ncol)

                # exponential function to fit, a*exp(-b*Av)
                def fexp(x, a, b):
                    return a * np.exp(-b * x)

                # fit exponential
                popt, pcov = curve_fit(fexp, avList, kphList)

                # store rate, 1/s
                self.phrates_extinction[db][k]["Draine1978"] = popt[1]

    # *******************
    # return rate (interface to dictionary), 1/s
    def getPhotoRate(self, database, process, radiation):
        return self.phrates[database][process][radiation]

    # ****************************
    # plot xsecs to PNG
    def plotXsec(self, pngFileName=None):

        import matplotlib.pyplot as plt

        # if name argument missing uses default
        if pngFileName is None:
            pngFileName = "pngs/xsec_" + self.nameFile + ".png"

        # turn off interactivity
        plt.ioff()

        # cancel current plot
        plt.clf()

        # set aesthetics
        plt.grid(b=True, color='0.65', linestyle='--')
        plt.xlabel("energy/eV")
        plt.ylabel("xsec/cm2")
        plt.title(self.name)

        hasPlot = False
        # loop on databases
        for (db, data) in self.xsecs.items():
            xdata = data["energy"]
            # loop on processes
            for (k, ydata) in data.items():
                # skip these keys
                if k in ["energy", "wavelength", "photoabsorption"]: continue
                # if no data no need to plot
                if len(ydata) > 0:
                    hasPlot = True
                    plt.plot(xdata, ydata, "-", label=k + " (" + db + ")")

        plt.legend(loc='best')
        if hasPlot: plt.savefig(pngFileName, dpi=150)

    # **********************
    def getHtmlName(self):
        name = list(self.name)
        if self.name.upper() == "E": return "e<sup>-</sup>"
        nameHtml = []
        for x in name:
            if utils.isNumber(x):
                y = "<sub>" + x + "</sub>"
            elif x == "+" or x == "-":
                y = "<sup>" + x + "</sup>"
            else:
                y = x
            nameHtml.append(y)
        return "".join(nameHtml)

    # ******************
    def getHrefName(self):
        url = "species_" + str(self.nameFile) + ".html"
        return "<a href=\"" + url + "\">" + self.getHtmlName() + "</a>"

    # **********************
    def getLatexName(self):
        name = list(self.name)
        if self.name.upper() == "E": return "e$^-$"
        nameLatex = []
        for x in name:
            if utils.isNumber(x):
                y = "$_" + x + "$"
            elif x == "+" or x == "-":
                y = "$^" + x + "$"
            else:
                y = x
            nameLatex.append(y)

        return "".join(nameLatex)

    # *****************
    # get polarizability, cm3
    def getPolarizability(self, polarizabilityData):

        # get polarizability from database
        if self.name in polarizabilityData:
            return polarizabilityData[self.name]
        return None

    # *****************
    # get "engineered" enthalpy kJ/mol
    def getEnthalpy(self, thermochemicalData, Tgas=298.15):

        # gas constant J/mol/K
        Rgas = 8.3144598

        # use so-called electron convention
        if self.name == "E": return 5. / 2. * Rgas * Tgas

        # CRs has no enthalpy
        if self.name == "CR": return 0e0

        # return None if unkonw element
        if self.name not in thermochemicalData:
            return None

        myData = thermochemicalData[self.name]
        if Tgas > myData["Tmid"]:
            coefs = myData["highT"]
        else:
            coefs = myData["lowT"]

        HRT = coefs[0] + \
              coefs[1] * Tgas / 2. + \
              coefs[2] * Tgas ** 2 / 3. + \
              coefs[3] * Tgas ** 3 / 4. + \
              coefs[4] * Tgas ** 4 / 5. + \
              coefs[5] / Tgas

        # kJ/mol
        return HRT * Rgas * Tgas / 1e3

    def getRadius(self, key):
        # return molecule radius in cm
        radii = {
            # Interatomic distance from Jeong et al 2000
            # DOI:10.1088/0953-4075/33/17/319
            "TiO2": 1.62e-8,
            # Interatomic distance O-Al-O (linear geometry) from
            # Archibong et al 1999 doi: 10.1021/jp983695n
            "Al2O3": 3.304e-8,
            # Half a bond length form Farrow et al 2014 doi:10.1039/C4CP01825G
            "MgO": 0.865e-8,
            # Half a bond length from Bromley et al 2016 doi:10.1039/c6cp03629e
            "SiO": 0.75765e-8
        }
        try:
            return radii[key]
        except KeyError:
            return None

    # *********************
    def makeHtmlPage(self, myNetwork):
        import urllib

        fname = "htmls/species_" + str(self.nameFile) + ".html"

        tableHeader = "<tr>" + ("<th>" * 30)

        # do xsec plot
        self.plotXsec()

        tableFormation = []
        tableDestruction = []
        for reactions in myNetwork.reactions:
            if self.name in [x.name for x in reactions.products]:
                tableFormation.append(reactions)
            if self.name in [x.name for x in reactions.reactants]:
                tableDestruction.append(reactions)

        fout = open(fname, "w")
        fout.write(utils.getFile("header.php"))
        fout.write("<p style=\"font-size:30px\">" + self.nameHtml + "</p>\n")
        fout.write("<br>\n")
        fout.write("<a href=\"indexSpecies.html\">back</a>\n")
        fout.write("<p>Enthalpy @ 298.15K: <b>" \
                   + str(self.getEnthalpy(myNetwork.thermochemicalData)) + "</b> kJ/mol</p>")
        polarizability = self.getPolarizability(myNetwork.polarizabilityData)
        if polarizability is not None: polarizability /= 1e-24  # cm3->AA3
        fout.write("<p>&alpha;: <b>" \
                   + str(polarizability) + "</b> &Aring;<sup>3</sup></p>")
        fnameURL = "species_allrates_" + str(self.nameFile) + ".html"
        fout.write("<p><a href=\"" + fnameURL + "\">All plots in a single page</a></p>")

        fout.write("<br><br>\n")
        fout.write("<p style=\"font-size:20px\">Formation channels</p><br>\n")
        fout.write("<table width=\"50%\">\n")
        fout.write(tableHeader + "\n")
        icount = 0
        for reaction in tableFormation:
            bgcolor = ""
            if icount % 2 != 0: bgcolor = utils.getHtmlProperty("tableRowBgcolor")
            fout.write(
                "<tr bgcolor=\"" + bgcolor + "\" valign=\"baseline\">" + reaction.getReactionHtmlRow(self) + "\n")
            icount += 1
        fout.write(tableHeader + "\n")
        fout.write("</table>\n")

        fout.write("<br><br>\n")
        fout.write("<p style=\"font-size:20px\">Destruction channels</p><br>\n")
        fout.write("<table width=\"50%\">\n")
        fout.write(tableHeader + "\n")
        icount = 0
        for reaction in tableDestruction:
            bgcolor = ""
            if icount % 2 != 0: bgcolor = utils.getHtmlProperty("tableRowBgcolor")
            fout.write(
                "<tr bgcolor=\"" + bgcolor + "\" valign=\"baseline\">" + reaction.getReactionHtmlRow(self) + "\n")
            icount += 1
        fout.write(tableHeader + "\n")
        fout.write("</table><br><br>\n")

        # put xsec png if data exists
        if len(self.xsecs) > 0:
            xsecPNG = "../pngs/xsec_" + self.nameFile + ".png"
            fout.write("<img src=\"" + xsecPNG + "\" width=\"500px\">\n")

            fname2 = "xsecs/phrateint_" + str(self.name) + ".txt"
            fout2 = open(fname2, 'w')  # for integrated photodestruction rates

            rows = []
            for db, data in self.phrates.items():
                # loop on rates
                for k, ydata in data.items():
                    for radName, kph in ydata.items():
                        bexp_str = ""
                        # include b from a*exp(-b*Av), WITHOUT SCATTERING!
                        if db in self.phrates_extinction:
                            if radName in self.phrates_extinction[db][k]:
                                bexp = self.phrates_extinction[db][k][radName]
                                bexp_str = str(round(bexp, 2))
                        tr = "<td>&nbsp;" + k + "<td>" + db + "<td>" + radName \
                             + "<td>" + utils.htmlExp(kph) + "<td>" \
                             + bexp_str + "\n"
                        rows.append([k + "_" + db + "_" + radName, tr])
                        line = k + " " + db + " " + radName + " " + str(kph) + " " + bexp_str
                        fout2.write(line + "\n")
            fout2.close()

            # table with integrated photorates
            thead = "<tr>" + "<th>" * 5
            fout.write("<br><br>\n")
            fout.write("<p style=\"font-size:20px\">Photochemical rates (1/s)</p><br>\n")
            fout.write("<table width=\"50%\">\n")
            fout.write(thead + "\n")
            btd = "b<sub>exp</sub><sup>*</sup>"
            fout.write("<tr><td>&nbsp;process<td>database<td>radiation<td>rate (1/s)<td>" + btd + "\n")
            fout.write(thead + "\n")
            icount = 0
            for row in sorted(rows, key=lambda x: x[0]):
                bgcolor = ""
                if icount % 2 != 0: bgcolor = utils.getHtmlProperty("tableRowBgcolor")
                fout.write("<tr bgcolor=\"" + bgcolor + "\">" + row[1])
                icount += 1
            fout.write(thead + "\n")
            fout.write("</table>\n")
            fout.write("<sup>*</sup>b from a*exp(-b*Av) extinction (WITHOUT SCATTERING!)\n")
            fout.write("<br><br>")

            link = "http://home.strw.leidenuniv.nl/~ewine/photo/index.php?file=display_species.php&species=" \
                   + urllib.parse.quote(self.name, safe='')
            fout.write("<a href=\"" + link + "\" target=\"_blank\">search on Leiden database</a><br>\n")
        else:
            link = "http://home.strw.leidenuniv.nl/~ewine/photo/data/photo_data/all_cross_sections/text_continuum/" \
                   + self.name + ".txt"
            fout.write("Cross-section missing: <a href=\"" + link \
                       + "\" target=\"_blank\">search on Leiden database</a> and copy to <i>xsecs/" \
                       + self.name + ".dat</i>\n")

        fout.write(utils.getFooter("footer.php"))
        fout.close()

    # ****************************
    def makeAllRatesHtmlPage(self, myNetwork, myOptions):

        fname = "htmls/species_allrates_" + str(self.nameFile) + ".html"

        tableFormation = []
        tableDestruction = []
        for reactions in myNetwork.reactions:
            if self.name in [x.name for x in reactions.products]:
                tableFormation.append(reactions)
            if self.name in [x.name for x in reactions.reactants]:
                tableDestruction.append(reactions)

        fout = open(fname, "w")
        fout.write(utils.getFile("header.php"))
        fout.write("<p style=\"font-size:30px\">" + self.nameHtml + "</p>\n")
        fout.write("<br>\n")
        fout.write("<a href=\"indexSpecies.html\">back</a>\n")

        reactions = {"formation": tableFormation, "destruction": tableDestruction}

        # loop on reactions block
        for k, data in reactions.items():
            fout.write("<br><br>\n")
            fout.write("<p style=\"font-size:20px\">" + k.title() + " channels</p><br>\n")
            # loop on reactions
            for reaction in data:
                # loop on variables
                for variable in myOptions.getRanges().keys():
                    fnamePNG = "../pngs/rate_" + str(reaction.getReactionHash()) + "_" + variable + ".png"
                    if reaction.hasVariable(myOptions, variable):
                        # fout.write("<img src=\""+fnamePNG+"\">")
                        linkURL = "<a href=\"rate_" + reaction.getReactionHash() + ".html\">details</a> for " + reaction.getVerbatimHtml()
                        fout.write(
                            "<img src=\"" + fnamePNG + "\" width=\"700px\" alt=\"&#9888; MISSING: " + reaction.getVerbatim() + "\"><br>" + linkURL + "<br><br>")

        fout.write(utils.getFooter("footer.php"))
        fout.close()
