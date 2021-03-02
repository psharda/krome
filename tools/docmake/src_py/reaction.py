import sys, species, utils, os, urllib
import ratefunctions
from copy import copy
from math import log10, log, exp, sqrt, pi


class reaction:
    index = -1
    verbatim = None
    reactants = []
    products = []
    Tmin = []
    Tmax = []
    rate = []
    verbatim = verbatimLatex = reactionHash = None
    nameHtml = verbatimHtml = reactionHtmlRow = None
    evaluation = []
    safeExtrapolate = dict()

    # ********************
    # parse csv reaction file row (constructor)
    def __init__(self, row, reactionFormat, atomSet, reactionType, speciesList, shortcuts=None, reference=None):

        if reactionType == "KIDA":
            self.parseFormatKIDA(row, reactionFormat, atomSet, reactionType, speciesList)
        elif reactionType == "UMIST":
            self.parseFormatUMIST(row, reactionFormat, atomSet, reactionType, speciesList)
        else:
            self.parseFormatKROME(row, reactionFormat, atomSet, reactionType, speciesList, shortcuts, reference)

    # ********************
    def parseFormatKROME(self, row, reactionFormat, atomSet, reactionType, speciesList, shortcuts=None, reference=None):

        if not reactionFormat.startswith("@format:"):
            sys.exit("ERROR: wrong format " + reactionFormat)
        splitFormat = reactionFormat.replace("@format:", "").split(",")
        splitFormat = [x.lower() for x in splitFormat]
        arow = row.strip().split(",", len(splitFormat) - 1)
        arow = [x.strip() for x in arow]

        specials = ["", "G"]

        self.rate = []
        self.reactants = []
        self.products = []
        self.Tmin = [None]
        self.Tmax = [None]
        self.reactionType = reactionType
        self.shortcuts = shortcuts
        self.reference = reference

        # loop on format parts (and parse species)
        for i in range(len(splitFormat)):
            part = splitFormat[i]
            # reaction index
            if part == "idx":
                self.index = int(arow[i])
            # reactant
            elif part == "r":
                if arow[i].upper() in specials: continue
                spec = self.linkOrCreateSpecies(arow[i], atomSet, speciesList)
                self.reactants.append(spec)
            # product
            elif part == "p":
                if arow[i].upper() in specials: continue
                spec = self.linkOrCreateSpecies(arow[i], atomSet, speciesList)
                self.products.append(spec)
            # min temperature
            elif part == "tmin":
                self.Tmin = [arow[i].replace("d", "e")]
                if arow[i].lower() == "none" or arow[i] == "": self.Tmin = [None]
            # max temperature
            elif part == "tmax":
                self.Tmax = [arow[i].replace("d", "e")]
                if arow[i].lower() == "none" or arow[i] == "": self.Tmax = [None]
            # rate coefficient
            elif part == "rate":
                self.rate.append(arow[i])
            else:
                print("ERROR: unknown format element " + part)
                sys.exit(reactionFormat)

        # add cosmic rays if not present
        hasCR = ("CR" in [x.name for x in self.reactants])
        if not hasCR and reactionType == "CR":
            spec = self.linkOrCreateSpecies("CR", atomSet, speciesList)
            self.reactants.append(spec)

        # check mass and charge conservation
        self.check()

    # ********************
    # replace names from KIDA style
    def speciesToKIDA(self, speciesName):

        if speciesName == "e-": speciesName = "E"
        if speciesName.startswith("p"): speciesName = speciesName[1:] + "_para"
        if speciesName.startswith("o"): speciesName = speciesName[1:] + "_ortho"
        if speciesName.startswith("m"): speciesName = speciesName[1:] + "_meta"

        return speciesName

    # ********************
    def parseFormatKIDA(self, row, reactionFormat, atomSet, reactionType, speciesList):

        specials = ["", "G", "PHOTON", "CR", "CRP"]

        # variables for cosmic rays and photochemistry
        CRvar = "user_crate"  # name of the CR flux variable
        Avvar = "user_Av"  # name of the Av variable

        # number of reactants and products expected in KIDA file
        maxReactants = 3
        maxProducts = 5

        # spacing format
        fmt = [11] * 3 + [1] + 5 * [11] + [1] + 3 * [11] + [8, 9] + [1, 4, 3] + 2 * [7] + [3, 6, 3, 2]
        # keys names
        keys = ["R" + str(i) for i in range(maxReactants)] + ["x"] \
               + ["P" + str(i) for i in range(maxProducts)] + ["x"] \
               + ["a", "b", "c"] + ["F", "g"] + ["x", "unc", "type"] \
               + ["tmin", "tmax"] + ["formula", "num", "subnum", "recom"]

        self.rate = []
        self.reactants = []
        self.products = []
        self.Tmin = [None]
        self.Tmax = [None]
        self.reactionType = reactionType

        srow = row.strip()

        dataRow = dict()
        position = 0
        # loop on format to get data from the row as a dictionary
        for i in range(len(fmt)):
            dataRow[keys[i]] = srow[position:position + fmt[i]].strip()
            startSpace = dataRow[keys[i]].startswith(" ")
            endSpace = dataRow[keys[i]].endswith(" ")
            hasSpace = (" " in dataRow[keys[i]])
            if not startSpace and not endSpace and hasSpace:
                print("ERROR: in KIDA network row element has spaces in the middle!")
                print(" Probably format problems:", dataRow[keys[i]])
                print(" Line here below")
                print(srow)
                sys.exit()
            position += fmt[i]

        self.Tmin = [dataRow["tmin"]]
        self.Tmax = [dataRow["tmax"]]

        for i in range(maxReactants):
            v = dataRow["R" + str(i)].strip()
            v = self.speciesToKIDA(v)
            if v.upper() in specials: continue
            spec = self.linkOrCreateSpecies(v, atomSet, speciesList)
            self.reactants.append(spec)

        for i in range(maxProducts):
            v = dataRow["P" + str(i)].strip()
            v = self.speciesToKIDA(v)
            if v.upper() in specials: continue
            spec = self.linkOrCreateSpecies(v, atomSet, speciesList)
            self.products.append(spec)

        # Formula is a number that referes to the formula needed to compute the rate coefficient of the reaction.
        # see http://kida.obs.u-bordeaux1.fr/help
        # 1: Cosmic-ray ionization (direct and undirect)
        # 2: Photo-dissociation (Draine)
        # 3: Kooij
        # 4: ionpol1
        # 5: ionpol2
        # 6: Troe fall-off (NOT SUPPORTED!)
        arow = dataRow
        arow["formula"] = int(arow["formula"])
        if arow["formula"] == 0:
            KK = "auto"
        elif arow["formula"] == 1:
            KK = arow["a"] + "*" + CRvar
        elif arow["formula"] == 2:
            KK = arow["a"]
            if float(arow["c"]) != 0e0: KK += "*exp(-" + arow["c"] + "*" + Avvar + ")"
        elif arow["formula"] == 3:
            KK = arow["a"]
            if float(arow["b"]) != 0e0: KK += "*(T32)**(" + arow["b"] + ")"
            if float(arow["c"]) != 0e0: KK += "*exp(-" + arow["c"] + "*invT)"
        elif arow["formula"] == 4:
            KK = arow["a"]
            if float(arow["b"]) != 1e0: KK += "*" + arow["b"]
            gpart = ""
            if float(arow["c"]) != 0e0: gpart = "+ 0.4767d0*" + arow["c"] + "*sqrt(3d2*invT)"
            KK += "*(0.62d0 " + gpart + ")"
        elif arow["formula"] == 5:
            KK = arow["a"]
            if float(arow["b"]) != 1e0: KK += "*" + arow["b"]
            gpart = ""
            if float(arow["c"]) != 0e0:
                gpart = "+ 0.0967d0*" + arow["c"] + "*sqrt(3d2*invT) + "
                gpart += arow["c"] + "**2*28.501d0*invT"
                KK += "*(1d0 " + gpart + ")"
        else:
            print("ERROR: KIDA formula " + str(arow["formula"]) + " not supported!")
            sys.exit()

        KK = KK.replace("--", "+").replace("++", "+").replace("-+", "-").replace("+-", "-")

        self.rate.append(KK)

    # ********************
    def parseFormatUMIST(self, row, reactionFormat, atomSet, reactionType, speciesList):

        # http://www.aanda.org/articles/aa/pdf/2013/02/aa20465-12.pdf
        # reaction no.:type:R1:R2:P1:P2:P3:P4:NE:[a:b:c:Tl:Tu:ST:ACC:REF]
        keys = "idx:type:R1:R2:P1:P2:P3:P4:NE:a:b:c:Tl:Tu:ST:ACC:REF".split(":")

        specials = ["", "CRP", "CRPHOT", "PHOTON"]

        # variables for cosmic rays and photochemistry
        CRvar = "user_crate"  # name of the CR flux variable
        Avvar = "user_Av"  # name of the Av variable
        refCRflux = 1.36e-17  # CR flux to scale variables

        self.rate = []
        self.reactants = []
        self.products = []
        self.Tmin = [None]
        self.Tmax = [None]
        self.reactionType = reactionType

        srow = row.strip()
        arow = [x.strip() for x in srow.split(":")]

        dataRow = dict()
        # loop on format to get data from the row as a dictionary
        for i in range(len(keys)):
            dataRow[keys[i]] = arow[i]

        self.Tmin = [dataRow["Tl"]]
        self.Tmax = [dataRow["Tu"]]

        for i in range(2):
            v = dataRow["R" + str(i + 1)].strip()
            if v == "e-": v = "E"
            if v.upper() in specials: continue
            spec = self.linkOrCreateSpecies(v, atomSet, speciesList)
            self.reactants.append(spec)
        for i in range(4):
            v = dataRow["P" + str(i + 1)].strip()
            if v == "e-": v = "E"
            if v.upper() in specials: continue
            spec = self.linkOrCreateSpecies(v, atomSet, speciesList)
            self.products.append(spec)

        if dataRow["type"] == "CP":
            KK = str(float(dataRow["a"]) / refCRflux) + "*" + CRvar
        elif dataRow["type"] == "CR":
            KK = str(float(dataRow["a"]) / refCRflux) + "*" + CRvar
            if float(dataRow["b"]) != 0e0: KK += "*(T32)**(" + dataRow["b"] + ")"
            KK += "*" + dataRow["c"] + "/(1d0-dustGrainAlbedo)"
            if float(dataRow["c"]) == 0e0: KK = "0d0"
        elif dataRow["type"] == "PH":
            KK = dataRow["a"]
            if float(dataRow["c"]) != 0e0: KK += "*exp(-" + dataRow["c"] + "*" + Avvar + ")"
        else:
            KK = dataRow["a"]
            if float(dataRow["b"]) != 0e0: KK += "*(T32)**(" + dataRow["b"] + ")"
            if float(dataRow["c"]) != 0e0: KK += "*exp(-" + dataRow["c"] + "*invT)"

        KK = KK.replace("--", "+").replace("++", "+").replace("-+", "-").replace("+-", "-")

        self.rate.append(KK)

    # **************
    # create a new species if not known, otherwise use known from list
    def linkOrCreateSpecies(self, speciesName, atomSet, speciesList):

        # create temporary species to compare parsed names
        tmpSpecies = species.species(speciesName, atomSet, skipPhotochemistry=True)

        # loop on species
        for spec in speciesList:
            # if name found return known
            if spec.name == tmpSpecies.name:
                return spec

        # create new species
        newSpecies = species.species(speciesName, atomSet)

        # append to list of species
        speciesList.append(newSpecies)

        return newSpecies

    # **************
    # check reaction charge and mass conservation
    def check(self):
        # check charge
        reactantsCharge = sum([x.charge for x in self.reactants])
        productsCharge = sum([x.charge for x in self.products])
        if reactantsCharge != productsCharge:
            print("ERROR: charge problems")
            print(self.getVerbatim())
            sys.exit()

        # check mass
        reactantsMass = sum([x.mass for x in self.reactants])
        productsMass = sum([x.mass for x in self.products])
        me = 9.10938356e-28  # electron mass, g
        # error if difference in mass is greater than half electron mass
        if abs(reactantsMass - productsMass) >= me / 2.:
            print("ERROR: mass problems")
            print(self.getVerbatim())
            print("Absolute difference (g):", abs(reactantsMass - productsMass))
            print("Mass reactants:")
            for RR in self.reactants:
                print(RR.name, RR.mass)
            print("Mass products:")
            for PP in self.products:
                print(PP.name, PP.mass)
            sys.exit()

    # ******************
    # update all reaction attributes
    def updateAttributes(self):
        self.getVerbatim()
        self.getVerbatimLatex()
        self.getVerbatimHtml()
        self.getReactionHash()

    # ********************
    # get verbatim reaction, e.g.H2+H2->H+H+H2
    def getVerbatim(self):
        if self.verbatim is not None: return self.verbatim
        reactantsName = [x.name for x in self.reactants]
        productsName = [x.name for x in self.products]
        self.verbatim = (" + ".join(reactantsName)) + " -> " + (" + ".join(productsName))
        return self.verbatim

    # ********************
    # get reaction verbatim in latex format, e.g. H$^+$+e$^-$$\to$H
    def getVerbatimLatex(self):
        if self.verbatimLatex is not None: return self.verbatimLatex
        reactantsName = [x.nameLatex for x in self.reactants]
        productsName = [x.nameLatex for x in self.products]
        self.verbatimLatex = (" + ".join(reactantsName)) + " $\\to$ " + (" + ".join(productsName))
        return self.verbatimLatex

    # ********************
    # get HTML verbatim, e.g. H<sub>2</sub>&rarr;H+H
    def getVerbatimHtml(self):
        if self.verbatimHtml is not None: return self.verbatimHtml
        reactantsName = [x.nameHtml for x in self.reactants]
        productsName = [x.nameHtml for x in self.products]
        self.verbatimHtml = (" + ".join(reactantsName)) + " &rarr; " + (" + ".join(productsName))
        return self.verbatimHtml

    # ********************
    def getRPHash(self, obj):
        rpName = sorted([x.nameFile for x in obj])
        return "_".join(rpName)

    # **************
    def getReactantsHash(self):
        return self.getRPHash(self.reactants)

    # **************
    def getProductsHash(self):
        return self.getRPHash(self.products)

    # ********************
    # get reaction unique hash, e.g. H_H__H2__standard
    def getReactionHash(self):
        if self.reactionHash is not None: return self.reactionHash
        self.reactionHash = self.getRPHash(self.reactants) + "__" + self.getRPHash(
            self.products) + "__" + self.reactionType
        return self.reactionHash

    # ********************
    # get reaction NON-unique (unsorted) hash, e.g. H_H__H2
    def getReactionHashUnsorted(self):
        reactantsName = [x.nameFile for x in self.reactants]
        productsName = [x.nameFile for x in self.products]
        return ("_".join(reactantsName)) + "__" + ("_".join(productsName))

    # ********************
    # get html table row with bold mySpecies when present
    def getReactionHtmlRow(self, mySpecies=None, mode=None):

        reactantsName = []
        # set bold to mySpecies name (reactants)
        for species in self.reactants:
            xspec = species.nameHref
            if mySpecies is not None:
                if species.name == mySpecies.name:
                    xspec = "<b>" + species.nameHtml + "</b>"
            reactantsName.append(xspec)

        productsName = []
        # set bold to mySpecies name (products)
        for species in self.products:
            xspec = species.nameHref
            if mySpecies is not None:
                if species.name == mySpecies.name:
                    xspec = "<b>" + species.nameHtml + "</b>"
            productsName.append(xspec)

        rpart = ("<td>+<td>".join(reactantsName))
        ppart = ("<td>+<td>".join(productsName))
        # fill empty spaces
        rspace = ("<td><td>".join([""] * (6 - len(reactantsName))))
        pspace = ("<td><td>".join([""] * (10 - len(productsName))))

        self.reactionHtmlRow = "<td>&nbsp;" + rpart + "<td>" + rspace + "<td>&rarr;<td>" + ppart + "<td>" + pspace
        self.reactionHtmlRow += "<td><a href=\"rate_" + self.getReactionHash() + ".html\">details</a>"

        # additional information depending on the required mode
        if mode == "joints":
            if len(self.evaluatedJoints) > 0:
                maxJointError = max([x["error"] for x in self.evaluatedJoints])
                warning = ("&#9888;" if (maxJointError > 1e-3) else "")
                self.reactionHtmlRow += "<td>" + str(round(maxJointError * 100, 2)) + "% " + warning
            else:
                self.reactionHtmlRow += "<td>N/A"
        return self.reactionHtmlRow

    # *******************
    def getAtoms(self):
        atoms = [x for x in self.getSpecies()]
        return list(set(atoms))

    # ********************
    # get list of species
    def getSpecies(self):
        return self.reactants + self.products

    # *******************
    # compute Langevin, cm3/s
    def computeLangevin(self, myNetwork):
        # Lanvevin only with two reactants
        if len(self.reactants) != 2: return None
        # needs a positive ion
        if +1 not in [x.charge for x in self.reactants]: return None
        # needs a neutral
        if 0 not in [x.charge for x in self.reactants]: return None

        polarizability = None
        # get neutral polarizability, cm3
        for reactant in self.reactants:
            if reactant.charge == 0:
                polarizability = reactant.getPolarizability(myNetwork.polarizabilityData)
        # need polarizability, cm3
        if polarizability is None:
            return None

        # inverse of reduced mass, 1/g
        invReducedMass = sum([1e0 / x.mass for x in self.reactants])

        elementaryCharge = 4.80320425e-10  # statC
        return 2e0 * pi * elementaryCharge * sqrt(polarizability * invReducedMass)

    # ********************
    # merge limits and rates with another rate
    def merge(self, myReaction):
        self.Tmin += myReaction.Tmin
        self.Tmax += myReaction.Tmax
        self.rate += myReaction.rate

    # *******************
    def evaluateRate(self, shortcuts, varRanges):
        self.evalRate(shortcuts, varRanges)
        self.evaluateExtrapolation(varRanges)
        self.evaluateJoints()
        self.saveEvals()

    # *******************
    # plot rate coefficient (alias)
    def plotRate(self, myOptions, pngFileName=None):
        self.doPlot(myOptions, pngFileName=pngFileName)

    # ********************
    # search for rate variables in the rate
    def getRateVariables(self, myOptions):
        rateVariables = []
        # loop on user defined variables
        for variable in myOptions.getRanges().keys():
            # loop on rates to search variable name
            for rate in self.rate:
                # append when variable found
                if self.hasSpecialRate:
                    rate = utils.getParentheticContents(rate, '()')[0][1]

                if variable.lower() in rate.lower():
                    rateVariables.append(variable)
        # if no variable found assumes Tgas
        if len(rateVariables) == 0:
            rateVariables = ["tgas"]

        # return unique list
        return list(set(rateVariables))

    # *****************
    # return rate reaction in KROME format
    def getKROMEformat(self, idx=1, includeTlimits=True):
        line = ""
        for i in range(len(self.rate)):
            fmt = "@format:idx,"
            fmt += ("R," * len(self.reactants))
            fmt += ("P," * len(self.products))
            if includeTlimits:
                if self.Tmin[i] is not None: fmt += "Tmin,"
                if self.Tmax[i] is not None: fmt += "Tmax,"
            fmt += "rate\n"
            RR = (",".join([x.name for x in self.reactants]))
            PP = (",".join([x.name for x in self.products]))
            Tmin = Tmax = ""
            if includeTlimits:
                if self.Tmin[i] is not None: Tmin = str(self.Tmin[i]) + ","
                if self.Tmax[i] is not None: Tmax = str(self.Tmax[i]) + ","
            data = str(idx + i + 1) + "," + RR + "," + PP + "," + Tmin + Tmax + self.rate[i] + "\n\n"
            line += fmt + data
        return line

    # *****************
    # boolean returned if variable found in rates
    def hasVariable(self, myOptions, variable):
        return variable.lower() in [x.lower() for x in self.getRateVariables(myOptions)]

    # ********************
    # evaluate rates
    def evalRate(self, shortcuts, varRanges):
        import re

        # number of points in the plot
        imax = 100

        self.evaluation = []
        self.evalRate = []
        self.warnings = []
        self.shortcutsFound = dict()
        self.rate2D = False

        # F90 expected operators in F90 rate expression
        ops = ["+", "-", "/", "*", "(", ")"]

        # loop on rates parts
        for icount in range(len(self.rate)):
            # get current rate
            rate = self.rate[icount]
            evaluation = dict()
            loopVariables = dict()  # variables present in rate
            vals = dict()  # variable range
            valsRange = dict()  # extra variable range for Tgas/extrapolation
            hasEval = dict()  # evaluation exist flag

            # check if rate is a function
            # that is defined in functionList
            # Currently only nucleation rate functions
            # TODO: extend with more functions
            # TODO: generalise for rates that are not only a function
            for specialRate in ratefunctions.functionList:
                if specialRate in rate:
                    self.hasSpecialRate = True
                    # name of the function
                    rateFunction = specialRate
                    # NOTE: currently it can only handle one function
                    break
                else:
                    self.hasSpecialRate = False

            if self.hasSpecialRate:
                # get list of the rate function arguments
                rateArguments = utils.getParentheticContents(rate, '()')[0][1].split(', ')

                for idx, arg in enumerate(rateArguments):
                    # replace "idx_species" with species object
                    if 'idx_' in arg:
                        speciesName = arg[4:]
                        for spec in self.getSpecies():
                            if spec.name == speciesName:
                                rateArguments[idx] = spec

                    # replace numbers with floats
                    elif utils.isNumber(arg):
                        rateArguments[idx] = utils.char2int(arg)

                    # surronds variables with "#"
                    else:
                        rateArguments[idx] = "#" + arg + "#"

                # remake rate string
                # objects also becomes a string (=unsuable)
                rate = re.sub(r'\(.*\)', '(' +
                              ", ".join([str(i) for i in rateArguments]) + ')', rate)

            else:
                substFound = True
                # loop until shortcut found and replaced
                while substFound:
                    # remove trailing comments (containing a "#")
                    rate = rate.split('#')[0]
                    # remove spaces
                    rate = rate.replace(" ", "").lower()

                    # surround F90 operators with # symbols
                    for op in ops:
                        rate = "#" + rate.replace(op, "#" + op + "#") + "#"

                    # remove double exponent operator
                    rate = rate.replace("#dexp#", "#exp#")

                    substFound = False
                    # split rate at #s
                    splitRate = [x for x in rate.split("#") if x != ""]
                    # sort shortcuts by size
                    var = sorted(shortcuts, key=lambda x: len(x), reverse=True)
                    # loop on shortcuts to replace
                    for variable in var:
                        # loop on rate parts
                        for i in range(len(splitRate)):
                            # if shortcut found replace with expression
                            if splitRate[i] == variable.lower():
                                splitRate[i] = "(" + shortcuts[variable] + ")"
                                self.shortcutsFound[variable] = shortcuts[variable]
                                substFound = True
                    # join rate back
                    rate = ("".join(splitRate))

                # replace F90 numbers for evaluation, d->e
                rate = rate.replace("d", "e")

                # surround F90 operators with # symbols
                for op in ops:
                    rate = "#" + rate.replace(op, "#" + op + "#") + "#"

            # check which/how many variables are in the rate
            # save the ones that are present
            for (variable, vrange) in varRanges.items():
                if "#" + variable.lower() + "#" in rate.lower():
                    loopVariables[variable] = vrange

            # if no variables present, skip rate
            if not loopVariables:
                continue
            # when only one variable, create a dummy for generalised 2D loop
            elif len(loopVariables) == 1:
                loopVariables['dummy'] = ['DUMMY']
            #
            elif len(loopVariables) == 2:
                self.rate2D = True

            else:
                print("Cannot handle %s number of variables in rate" % (len(loopVariables)))

            # loop on available ranges
            for variable, vrange in loopVariables.items():
                hasEval[variable] = False
                # check if Tgas
                isTgas = (variable.lower() == "tgas")
                if variable == 'dummy':
                    vals[variable] = ['']
                    continue
                # get range limits
                (varMin, varMax) = vrange
                # log limits
                logVarMin = log10(varMin)
                logVarMax = log10(varMax)
                # create variable range
                vals[variable] = [i * (logVarMax - logVarMin) / (imax - 1) + logVarMin for i in range(imax)]
                vals[variable] = [1e1 ** x for x in vals[variable]]

                # when Tgas check add limits if any
                if isTgas:
                    Tmin = varMin
                    Tmax = varMax
                    if self.Tmin[icount] is not None: Tmin = float(utils.replaceTlims(self.Tmin[icount]))
                    if self.Tmax[icount] is not None: Tmax = float(utils.replaceTlims(self.Tmax[icount]))
                    Tmin = max(varMin, Tmin)
                    Tmax = min(varMax, Tmax)
                    valsRange[variable] = [x for x in vals[variable] if (Tmin <= x and x <= Tmax)]
                    valsRange[variable] = [Tmin] + valsRange[variable] + [Tmax]

                    # add additional points close to the limits (needed by evaluate joints)
                    vals[variable] += [Tmin, Tmax]
                    for dx in [0.5, 1e0]:
                        if Tmin - dx > 0: vals[variable] += [Tmin - dx]
                        vals[variable] += [Tmin + dx, Tmax - dx, Tmax + dx]

                    # sort Tgas values
                    vals[variable] = sorted(vals[variable])

                # store evaluated rate
                self.evalRate.append(rate)
                # store xdata and init ydata
                evaluation[variable] = dict()
                evaluation[variable]["xdata"] = vals[variable]
                evaluation[variable]["ydata"] = []
                # if Tgas store as interval
                if isTgas:
                    evaluation[variable]["xdataRange"] = valsRange[variable]
                    evaluation[variable]["ydataRange"] = []

                if self.rate2D:
                    evaluation[variable]["zdata"] = []
                    if isTgas:
                        evaluation[variable]["zdataRange"] = []

            # evaluate rate full range
            # loop over ranges of all found variables
            keyVars = list(loopVariables.keys())
            for valFirst in vals[keyVars[0]]:
                for valSecond in vals[keyVars[1]]:
                    # check if rate is a function
                    if self.hasSpecialRate:
                        # copy list to avoid replacing arguments permanetly
                        rateArgumentsNew = rateArguments[:]
                        for idx, arg in enumerate(rateArguments):
                            # replace variable with its float value
                            try:
                                arg = arg.replace("#" + keyVars[0] + "#", str(valFirst))
                                rateArgumentsNew[idx] = float(arg.replace("#" + keyVars[1] + "#", str(valSecond)))
                            except:
                                continue
                        # call the special rate function
                        # NOTE: I tried using 'exec(rate)'
                        # but python finds this illegal
                        # so an ugly if elif for all functions
                        # in the solution...
                        if rateFunction == "cluster_growth_rate":
                            yvalue = ratefunctions.cluster_growth_rate(rateArgumentsNew[0],
                                                                       rateArgumentsNew[1],
                                                                       rateArgumentsNew[2])
                            # save evaluation flag
                            # I know it looks ugly... :/ => TODO
                            if self.rate2D:
                                hasEval[keyVars[0]] = True
                                hasEval[keyVars[1]] = True
                            elif keyVars[0] == 'dummy':
                                hasEval[keyVars[1]] = True
                            else:
                                hasEval[keyVars[0]] = True
                        elif rateFunction == "cluster_destruction_rate":
                            yvalue = ratefunctions.cluster_destruction_rate(rateArgumentsNew[0],
                                                                            rateArgumentsNew[1],
                                                                            rateArgumentsNew[2])
                            # save evaluation flag
                            if self.rate2D:
                                hasEval[keyVars[0]] = True
                                hasEval[keyVars[1]] = True
                            elif keyVars[0] == 'dummy':
                                hasEval[keyVars[1]] = True
                            else:
                                hasEval[keyVars[0]] = True
                        elif rateFunction == "steady_state_nucleation_rate":
                            yvalue = ratefunctions.steady_state_nucleation_rate(rateArgumentsNew[0],
                                                                                rateArgumentsNew[1],
                                                                                rateArgumentsNew[2],
                                                                                rateArgumentsNew[3])
                            # save evaluation flag
                            if self.rate2D:
                                hasEval[keyVars[0]] = True
                                hasEval[keyVars[1]] = True
                            elif keyVars[0] == 'dummy':
                                hasEval[keyVars[1]] = True
                            else:
                                hasEval[keyVars[0]] = True

                        else:
                            yvalue = None
                            print('ERROR: %s not defined in rateFunction.py' % rateFunction)

                    else:
                        # replace variables with their floats
                        k = rate.replace("#" + keyVars[0].lower() + "#", str(valFirst))
                        k = k.replace("#" + keyVars[1].lower() + "#", str(valSecond)).replace("#", "")
                        try:
                            yvalue = eval(k)
                            if self.rate2D:
                                hasEval[keyVars[0]] = True
                                hasEval[keyVars[1]] = True
                            elif keyVars[0] == 'dummy':
                                hasEval[keyVars[1]] = True
                            else:
                                hasEval[keyVars[0]] = True
                        except:
                            yvalue = None
                    if self.rate2D:
                        # save in one of two variable dicts
                        evaluation[keyVars[0]]["zdata"].append(yvalue)
                    else:
                        # save in variable dicts that is not 'dummy'
                        if keyVars[0] == 'dummy':
                            evaluation[keyVars[1]]["ydata"].append(yvalue)
                        else:
                            evaluation[keyVars[0]]["ydata"].append(yvalue)

                # #check if negative values found
                # if not self.rate2D:
                #     if( min(evaluation[variable]["ydata"]) < 0 and hasEval[variable]):
                #         if(isTgas):
                #             self.warnings.append("negative values when extrapolated")
                #         else:
                #             self.warnings.append("negative values found")

            # no extrapolation for 2D color plots
            if not self.rate2D:
                for variable, vrange in loopVariables.items():
                    # check if Tgas
                    isTgas = (variable.lower() == "tgas")
                    # evaluate rate limited range
                    if isTgas and hasEval[variable]:
                        if self.hasSpecialRate:
                            # replace variable with its float value
                            for idx, arg in enumerate(rateArguments):
                                try:
                                    argmin = float(arg.replace("#" + variable + "#", str(Tmin)))
                                    argmax = float(arg.replace("#" + variable + "#", str(Tmax)))
                                except:
                                    continue

                            # call special rate function
                            # same note as above
                            if rateFunction == "cluster_growth_rate":
                                kmin = ratefunctions.cluster_growth_rate(rateArguments[0],
                                                                         rateArguments[1],
                                                                         argmin)
                                kmax = ratefunctions.cluster_growth_rate(rateArguments[0],
                                                                         rateArguments[1],
                                                                         argmax)

                            elif rateFunction == "cluster_destruction_rate":
                                kmin = ratefunctions.cluster_destruction_rate(rateArguments[0],
                                                                              rateArguments[1],
                                                                              argmin)

                                kmax = ratefunctions.cluster_destruction_rate(rateArguments[0],
                                                                              rateArguments[1],
                                                                              argmax)
                            else:
                                print('%s not defined in rateFunction.py' % rateFunction)

                        else:
                            try:
                                kmin = eval(rate.replace("#" + variable.lower() + "#", str(Tmin)).replace("#", ""))
                                kmax = eval(rate.replace("#" + variable.lower() + "#", str(Tmax)).replace("#", ""))
                            except:
                                print("ERROR: problem evaluating rate at limits")
                                print("limits:", str(Tmin), str(Tmax))
                                print("rate:", rate.replace("#", ""))
                                sys.exit()

                        evaluation[variable]["xlimits"] = [Tmin, Tmax]
                        evaluation[variable]["ylimits"] = [kmin, kmax]
                        for val in valsRange[variable]:
                            if self.hasSpecialRate:
                                # copy list t avoid permanent replacements
                                rateArgumentsNew = rateArguments[:]
                                # replace variable with its float value
                                for idx, arg in enumerate(rateArguments):
                                    try:
                                        rateArgumentsNew[idx] = float(arg.replace("#" + variable + "#", str(val)))
                                    except:
                                        continue

                                # call special rate function
                                # same note as above
                                if rateFunction == "cluster_growth_rate":
                                    yvalue = ratefunctions.cluster_growth_rate(rateArgumentsNew[0],
                                                                               rateArgumentsNew[1],
                                                                               rateArgumentsNew[2])
                                    hasEval[variable] = True

                                elif rateFunction == "cluster_destruction_rate":
                                    yvalue = ratefunctions.cluster_destruction_rate(rateArgumentsNew[0],
                                                                                    rateArgumentsNew[1],
                                                                                    rateArgumentsNew[2])
                                    hasEval[variable] = True
                                else:
                                    yvalue = None
                                    print('%s not defined in rateFunction.py' % rateFunction)

                            else:
                                k = rate.replace("#" + variable.lower() + "#", str(val)).replace("#", "")
                                try:
                                    yvalue = eval(k)
                                    hasEval[variable] = True
                                except:
                                    yvalue = None

                            evaluation[variable]["ydataRange"].append(yvalue)

                        # check if negative values found (when limited range)
                        if min(evaluation[variable]["ydataRange"]) < 0 and hasEval[variable]:
                            self.warnings.append("negative values found")
            for variable, vrange in loopVariables.items():
                # if not evaluated put none
                if not hasEval[variable]:
                    # add warning if rate not evaluated
                    self.warnings.append("no rate evaluation")
                    evaluation[variable] = None

            # store as class attribute
            self.evaluation.append(evaluation)

    # **************************
    # do plot (PNG)
    def doPlot(self, myOptions, pngFileName=None):
        import matplotlib
        import numpy as np
        import matplotlib.colors as colors
        # try to load AGG for PNG rendering (slightly faster)
        # try:
        #    matplotlib.use('AGG')
        # except:
        #    pass

        import matplotlib.pyplot as plt

        # turn off interactivity
        plt.ioff()

        # cancel current plot
        plt.clf()
        # max orders of magnitude y axis
        yspanMax = 1e-20
        hasPlot = False
        ydataAll = []
        # variable names to save plots for
        saveVariables = []

        # check if rate has two variables
        if self.rate2D:
            # loop on range varibles
            for rng in myOptions.range:
                # get range name
                variable = rng.split("=")[0].strip()
                # loop on different limited ranges
                for evaluation in self.evaluation:
                    # skip variables that are not in the rate
                    if variable not in evaluation:
                        continue
                    data = evaluation[variable]

                    if data['zdata']:
                        zdata = data['zdata']
                        ydata = data['xdata']
                        plt.ylabel(variable)

                    else:
                        xdata = data['xdata']
                        plt.xlabel(variable)

                    hasPlot = True
                    saveVariables.append(variable)

            # adapt data to enable 2D color plot
            Nrow = len(ydata)
            Ncol = len(xdata)
            zdata = np.reshape(zdata, (Nrow, Ncol))

            # set ranges
            zmin = max(zdata.max() * yspanMax, zdata.min())
            zmax = zdata.max()
            xRange = max(xdata) / min(xdata)
            yRange = max(ydata) / min(ydata)
            zRange = zmax / zmin

            # with data range is too large, make logaritmic color plot
            if zRange > 10:
                plt.pcolormesh(xdata, ydata, zdata, cmap='viridis', rasterized=True,
                               norm=colors.LogNorm(vmin=zmin, vmax=zdata.max()))
            else:
                plt.pcolormesh(xdata, ydata, zdata, cmap='viridis', rasterized=True)

            plt.colorbar(label='Reaction rate', extend='min')

            # make logaritmic when range is too large
            if xRange > 99:
                plt.xscale('log')
            if yRange > 99:
                plt.yscale('log')

            plt.title(self.getVerbatimLatex())

            # if argument is not present automatic file
            if pngFileName is None and hasPlot:
                # make two plots for 2D rats so they appear twice in the html
                for var in saveVariables:
                    pngFileName = "pngs/rate_" + str(self.getReactionHash()) + "_" + var + ".png"

            # plot only if data are available
            if hasPlot and not os.path.exists(pngFileName):
                # if value found save plot to png file
                plt.savefig(pngFileName, dpi=150)

        else:

            # loop on range varibles
            for rng in myOptions.range:
                # get range name
                variable = rng.split("=")[0].strip()
                plt.clf()
                # loop on different limited ranges
                for evaluation in self.evaluation:
                    if (not (variable in evaluation)): continue
                    data = evaluation[variable]
                    if (data == None): continue
                    xdata = data["xdata"]
                    ydata = data["ydata"]
                    if all([yd == 0.0 for yd in ydata]):
                        print("WARNING: The rate for {0} is zero; skipping the plot!".format(self.getVerbatim()))
                        continue  # all rate data are zero, so skip plotting
                    hasPlot = True

                    if variable.lower() == "tgas":
                        # plot full range
                        plt.loglog(xdata, ydata, "r--")

                        # if Tgas use limited ranges and plot limit points
                        xdataRange = data["xdataRange"]
                        ydataRange = data["ydataRange"]
                        ydataAll += ydataRange
                        plt.loglog(evaluation[variable]["xlimits"], evaluation[variable]["ylimits"], "ro")
                        plt.loglog(xdataRange, ydataRange, "b-")
                    else:
                        xdataRange = xdata
                        ydataRange = ydata
                        ydataAll += ydataRange
                        plt.loglog(xdataRange, ydataRange)

                # if argument is not present automatic file
                if pngFileName is None and hasPlot:
                    pngFileName = "pngs/rate_" + str(self.getReactionHash()) + "_" + variable + ".png"

                # plot only if data are available
                if hasPlot and not os.path.exists(pngFileName):
                    plt.grid(b=True, color='0.65', linestyle='--')
                    # plot limited range
                    plt.xlabel(variable)
                    plt.ylabel("rate")
                    plt.title(self.getVerbatimLatex())
                    # set limits including max span
                    plt.ylim(max(max(ydataAll) * yspanMax, min(ydataAll) * 1e-1), max(ydataAll) * 1e1)
                    # set limits if constant
                    if min(ydataAll) == max(ydataAll): plt.ylim(max(ydataAll) * 1e-1, max(ydataAll) * 1e1)

                    # if value found save plot to png file
                    plt.savefig(pngFileName, dpi=150)

    # ******************
    # evaluate rate extrapolation for the current reaction
    def evaluateExtrapolation(self, varRanges):

        # skip when rate has two variables
        if self.rate2D:
            return

        # init Tgas limits
        xMin = 1e99
        xMax = -1e99
        # init flags
        hasData = isIncreasingMin = isDecreasingMax = False
        isAlwaysPositiveMin = isAlwaysPositiveMax = False
        # loop on different limited ranges
        for evaluation in self.evaluation:
            data = dict()
            # get only Tgas data
            for variable, vdata in evaluation.items():
                if variable.lower() != "tgas":
                    continue
                # store data
                data = vdata
                # store ranges
                varRange = varRanges[variable]

            # skip missing data
            if not data: continue
            hasData = True
            # copy data locally (evaluation in the rate Tgas range)
            xdataRange = data["xdataRange"]
            ydataRange = data["ydataRange"]
            # copy data locally (evaluation in the whole Tgas range)
            xdata = data["xdata"]
            ydata = data["ydata"]
            # number of data points
            ndata = len(xdata)

            # check smaller rate interval
            if min(xdataRange) < xMin:
                # store min value
                xMin = min(xdataRange)
                # get ydata outside interval
                ydataOutside = [ydata[i] for i in range(ndata) if (xdata[i] < xMin)]
                # check if ydata are increasing outside
                isIncreasingMin = (sorted(ydataOutside) == ydataOutside)
                # check if data are always positive outside (only if data are present)
                isAlwaysPositiveMin = True
                if len(ydataOutside) > 0: isAlwaysPositiveMin = (min(ydataOutside) > 0e0)
                # store min Tgas in data structure
                self.safeExtrapolate["Tmin"] = xMin

            if max(xdataRange) > xMax:
                xMax = max(xdataRange)
                ydataOutside = [ydata[i] for i in range(ndata) if (xdata[i] > xMax) and ydata[i] is not None]
                isDecreasingMax = (sorted(ydataOutside) == ydataOutside[::-1])
                isAlwaysPositiveMax = True
                if len(ydataOutside) > 0:
                    isAlwaysPositiveMax = min(ydataOutside) > 0e0
                self.safeExtrapolate["Tmax"] = xMax

        # check extrapolation only if has data
        if hasData:
            # store extrapolated Tgas limits
            self.safeExtrapolate["TminExtrapolated"] = min(varRange)
            self.safeExtrapolate["TmaxExtrapolated"] = max(varRange)
            # store if extrapolation is safe or not
            self.safeExtrapolate["lower"] = (isAlwaysPositiveMin and isIncreasingMin)
            self.safeExtrapolate["upper"] = (isAlwaysPositiveMax and isDecreasingMax)

    # ****************
    # evaluate rate joints
    def evaluateJoints(self):

        self.evaluatedJoints = []

        dataAll = []
        # loop on different limited ranges
        for evaluation in self.evaluation:
            # get only Tgas data
            for (variable, vdata) in evaluation.items():
                if variable.lower() != "tgas": continue
                if vdata is None: continue
                # store data
                dataAll.append(vdata)

        # if less than two intervals ignore
        if len(dataAll) < 2: return

        # min distance to determine close limits
        distanceThreshold = 2e0
        # temperature shift
        dx = 1e0
        # loop on ranges
        for data1 in dataAll:
            # loop on ranges
            for data2 in dataAll:
                # check distance
                if abs(data1["xdataRange"][-1] - data2["xdataRange"][0]) < distanceThreshold:
                    # try to pick the limit in the extrapolated other range
                    # otherwise shift by 1K
                    try:
                        idx2 = data2["xdata"].index(data1["xdataRange"][-1])
                    except:
                        idx2 = data2["xdata"].index(data1["xdataRange"][-1] + dx)

                    # store data
                    yrate = data1["ydataRange"][-1]
                    yeval = data2["ydata"][idx2]
                    xeval = data2["xdata"][idx2]
                    joint = dict()
                    joint["limit1"] = [data1["xdataRange"][-1], data1["ydataRange"][-1]]
                    joint["limit2"] = [data2["xdataRange"][0], data2["ydataRange"][0]]
                    joint["extrapolation"] = [xeval, yeval]
                    joint["error"] = abs(yrate - yeval) / (yrate + 1.0e-99)
                    self.evaluatedJoints.append(joint)

    # ****************
    # save evaluation as a json structure
    def saveEvals(self):
        import json

        # json file name
        fname = "evals/rate_" + str(self.getReactionHash()) + ".json"

        # convert to json
        jsonData = json.dumps(self.evaluation)

        # dump to file
        fout = open(fname, "w")
        fout.write(jsonData)
        fout.close()

    # ****************
    # make corresponding HTML page
    def makeHtmlPage(self, myOptions, myNetwork):

        fname = "htmls/rate_" + str(self.getReactionHash()) + ".html"

        reactantsList = [x.getHtmlName() for x in self.reactants]
        productsList = [x.getHtmlName() for x in self.products]

        header = "<tr><th><th><th>\n"

        table = []
        tableNotes = []
        # table.append(["reactants", (", ".join(reactantsList))])
        # table.append(["products", (", ".join(productsList))])
        for icount in range(len(self.rate)):
            table.append(["header", header])
            table.append(["rate", self.rate[icount]])
            table.append(["Tmin", self.Tmin[icount]])
            table.append(["Tmax", self.Tmax[icount]])

        for (shortcutName, shortcutExpression) in self.shortcutsFound.items():
            tableNotes.append(shortcutName + " = " + shortcutExpression)

        for warning in self.warnings:
            tableNotes.append(warning)

        fout = open(fname, "w")
        fout.write(utils.getFile("header.php"))
        fout.write("<p style=\"font-size:30px\">" + self.getVerbatimHtml() + "</p>\n")
        fout.write("<br>\n")
        fout.write("<a href=\"indexReactions.html\">back</a><br>\n")
        urlencoded = urllib.parse.quote_plus(" + ".join([x.name for x in self.reactants]))
        urlkida = "http://kida.obs.u-bordeaux1.fr/search.html?species=" + urlencoded \
                  + "&reactprod=reactants&astroplaneto=Both&ionneutral=ion&isomers=1&ids="
        fout.write("<a href=\"" + urlkida + "\" target=\"_blank\">search in KIDA</a><br>\n")
        urlJSON = "../evals/rate_" + str(self.getReactionHash()) + ".json"
        fout.write("<a href=\"" + urlJSON + "\">get rate evaluation in JSON format</a>\n")
        fout.write("<br><br>\n")

        langevinRate = self.computeLangevin(myNetwork)
        if langevinRate is not None: fout.write("Langevin rate: " + utils.htmlExp(langevinRate) \
                                                + " cm<sup>3</sup>s<sup>-1</sup><br>\n")

        allSpecies = sorted(self.getSpecies(), key=lambda x: x.name)
        hrefNames = [x.getHrefName() for x in allSpecies]
        fout.write("Species involved: " + (", ".join(hrefNames)) + "<br><br>")

        fout.write("<table>\n")
        for (label, value) in table:
            if value is None: continue
            if label == "header":
                fout.write(value)
            else:
                separator = "&nbsp;&nbsp;:&nbsp;&nbsp;"
                fout.write("<tr><td>" + label + "<td>" + separator + "<td>" + value + "\n")
        fout.write(header)
        fout.write("</table><br><br>\n")

        bulletPoint = "&nbsp;&nbsp;&#9656;&nbsp;"
        # add notes and warnings
        if len(tableNotes) > 0:
            fout.write("Notes and warnings:<br>\n")
            for value in tableNotes:
                if value is not None: fout.write(bulletPoint + value + "<br>\n")

        # extrapolation lower limit
        if "TminExtrapolated" in self.safeExtrapolate:
            TminExtrapolated = utils.htmlExpBig(self.safeExtrapolate["TminExtrapolated"])
            Tmin = utils.htmlExpBig(self.safeExtrapolate["Tmin"])
            extrapolCheckMin = ("SAFE" if self.safeExtrapolate["lower"] else "NOT SAFE")
            if self.safeExtrapolate["TminExtrapolated"] != self.safeExtrapolate["Tmin"]:
                fout.write(bulletPoint + "Extrapolation in range [" + TminExtrapolated \
                           + ", " + Tmin + "] K is <b>" + extrapolCheckMin + "</b><br>")

        # extrapolation upper limit
        if "TmaxExtrapolated" in self.safeExtrapolate:
            TmaxExtrapolated = utils.htmlExpBig(self.safeExtrapolate["TmaxExtrapolated"])
            Tmax = utils.htmlExpBig(self.safeExtrapolate["Tmax"])
            extrapolCheckMax = ("SAFE" if self.safeExtrapolate["upper"] else "NOT SAFE")
            if self.safeExtrapolate["TmaxExtrapolated"] != self.safeExtrapolate["Tmax"]:
                fout.write(bulletPoint + "Extrapolation in range [" + Tmax + ", " \
                           + TmaxExtrapolated + "] K is <b>" + extrapolCheckMax + "</b><br>")

        # JOINTS evaluation
        if len(self.evaluatedJoints) > 0:
            sep = "<td>&nbsp;:&nbsp;<td>"
            header = "<tr><th><th><th><th>"
            fout.write("<br>")
            fout.write("Evaluated joints:<br>")
            fout.write("<table>")
            fout.write(header)
            fout.write("<tr><td><td><td>Tgas/K<td>rate")
            for joint in self.evaluatedJoints:
                fout.write(header)
                fout.write("<tr><td>limit1 (Tgas,rate)" + sep \
                           + ("<td>".join([str(x) for x in joint["limit1"]])))
                fout.write("<tr><td>limit2 (Tgas,rate)" + sep \
                           + ("<td>".join([str(x) for x in joint["limit2"]])))
                fout.write("<tr><td>extrapolated (Tgas,rate)" + sep \
                           + ("<td>".join([str(x) for x in joint["extrapolation"]])))
                warning = ""
                if joint["error"] > 1e-3: warning = "&#9888;"
                fout.write("<tr><td>error" + sep + str(round(joint["error"] * 100, 2)) + "% " + warning)
            fout.write(header)
            fout.write("</table>")

        # PLOT
        for rng in myOptions.range:
            (rangeName, rangeValue) = [x.strip() for x in rng.split("=")]
            plotFileName = "pngs/rate_" + str(self.getReactionHash()) + "_" + rangeName + ".png"

            hasPlot = False
            # loop on different limits to check if data are present
            for evaluation in self.evaluation:
                if rangeName not in evaluation: continue
                data = evaluation[rangeName]
                if data is None: continue
                if self.rate2D:
                    ratedata = data["zdata"]
                else:
                    ratedata = data["ydata"]
                if all([yd == 0.0 for yd in ratedata]):
                    fout.write(bulletPoint + "The rate for this reaction is <b>ZERO</b><br>")
                    continue
                if evaluation[rangeName] is not None:
                    hasPlot = True
                    break

            if hasPlot:
                fout.write("<img width=\"700px\" src=\"../" + plotFileName + "\">\n")

        fout.write(utils.getFooter("footer.php"))
        fout.close()

    # ********************
    # make a LaTeX format of reaction
    # TODO: remove deprecated arguments
    def reaction2latex(self, temperatureShortcuts, variableShortcuts, deferredShortcuts,
                       cntMergedReactions, idxMerged, cntTotalReactions, cntAllReactions):
        import sys
        import re
        from options import latexoptions as opts
        # latex format uses \usepackage{chemformula} in LaTeX
        # e.g. \ch{H2 + H -> H + H + H}

        # index of unique reactions
        idxUniqueReaction = str(cntTotalReactions - cntMergedReactions)

        if idxMerged == 0:
            # LaTeX index
            if opts.sorted_by == 'alphabetic':
                idxTex = idxUniqueReaction
            else:
                idxTex = str(self.uid)

            # latex reaction
            reactionTex = "\\ch{"
            # loop over reactants
            for r in self.reactants:
                spec = r.name
                if spec == "E":
                    spec = "e-"
                # remove CR as species in LaTeX format
                # this can be changed if the user prefers otherwise
                if spec == "CR":
                    continue

                # chemformula requires species with more than one ionizations to have them bracketed
                if "++" in spec:
                    spec = re.sub(r"(\+{2,})", r"^{\1}", spec)
                if "--" in spec:
                    spec = re.sub(r"(\-{2,})", r"^{\1}", spec)

                reactionTex += spec + " + "
            # if photo-reaction, add photon as reactant
            if self.reactionType == "photo":
                reactionTex += "$\\gamma$"
            else:
                # remove trailing " + "
                reactionTex = reactionTex[:-3]
            # if cosmic ray reaction, make LaTeX format
            # e.g. \ch{H2 ->[CR] H + H}
            if self.reactionType == "CR":
                reactionTex += " ->[CR] "
            elif self.reactionType == "catalysis":
                reactionTex += " + M -> "
            else:
                reactionTex += " -> "

            # loop over products
            for p in self.products:
                spec = p.name
                if spec == "E":
                    spec = "e-"
                if "++" in spec:
                    spec = re.sub(r"(\+{2,})", r"^{\1}", spec)
                if "--" in spec:
                    spec = re.sub(r"(\-{2,})", r"^{\1}", spec)
                reactionTex += spec + " + "

            # add generic catalysist to catalysis reaction
            if self.reactionType == "catalysis":
                reactionTex += "M"
            # add photon to reaction with only one product
            elif len(self.products) == 1:
                reactionTex += "$\\gamma$"
            else:
                # remove trailing " + "
                reactionTex = reactionTex[:-3]

            reactionTex += "}"

            # default LaTeX reference as a LaTeX command
            if self.reference is None:
                refTex = "\\dbDefault{}"
            else:
                refTex = str(self.reference[0])

        else:
            idxTex = ""
            reactionTex = ""
            refTex = ""

        # look in database if this is an auto rate
        if self.rate[idxMerged] == "auto":
            self.rate[idxMerged], Emin, Emax = self.find_auto_rate(self.reactants, self.products, self.reactionType)

        # LaTeX rate
        rateTex, message, numlines = self.rate2latex(self.rate[idxMerged], temperatureShortcuts, variableShortcuts,
                                                     deferredShortcuts)
        rateTex = "$" + rateTex + "$"

        # LaTeX temperature limits
        limitsTex = self.tempRange2latex(idxMerged)
        limitsTex = "$" + limitsTex + " $"

        # k symbol
        if idxMerged == 0:
            kTex = "k$_{" + idxTex + "}$"
        else:
            kTex = ""

        return [idxTex, reactionTex, kTex, rateTex, limitsTex, refTex], message, numlines

    # ********************
    # make a LaTeX format of reaction
    def rate2latex(self, rate, temperatureShortcuts, variableShortcuts, deferredShortcuts):
        import re
        from options import latexoptions as opts
        if opts.latex_backend == "pytexit":
            import pytexit
        elif opts.latex_backend == "sympy":
            import sympy as sp
            num = sp.__version__.count('.') - 1
            sp_version = float(sp.__version__.rsplit('.', num)[0])
            if sp_version >= 1.3:
                print("ERROR: The LaTeX conversion currently only works with and older"
                      " version of SymPy (<1.3). Symbols no longer automatically"
                      " convert to functions when called."
                      )
                sys.exit()
        else:
            print("WARNING: Option '{}' is not reconized as LaTeX convertor"
                  + "Adapt the '{}' file").format(opts.latex_backend, opts.__name__)
            exit()

        debug = False
        maxRateLength = 100
        message = ""  # optional warning

        symboltable = utils.getSymbolTable()

        originalRate = copy(rate)

        # put all variables with corresponding values in rate
        if variableShortcuts:
            rate = utils.replaceShortcuts(rate, variableShortcuts, deferredShortcuts.keys())

        # replace shortcuts, loop needs to be reversed order for variable dependencies
        # skip T32 and Te to keep as symbol
        rate = utils.replaceShortcuts(rate, temperatureShortcuts[1:], deferredShortcuts.keys())

        # make sympy friendly
        rate = utils.replaceFortranVar("dexp", "exp", rate)
        rate = utils.replaceFortranVar("log", "ln", rate)
        rate = utils.replaceFortranVar("log10", "log", rate)
        # Double prec exp notation to use e, e.g. 2d3 -> 2e3
        rate = re.sub("([0-9]*\.*[0-9]*)d([+-]*[0-9]{1,})", r"\1e\2", rate)
        # Remove double prec suffix, e.g. 1.0_8 -> 1.0
        rate = re.sub("([0-9])_8", r"\1", rate)
        rate = re.sub("\\*1e0/([a-zA-Z][a-zA-Z0-9_]*)", "/\\1",
                      rate)  # Replace *1/x by /x to avoid dangling *1 in fractions
        rate = re.sub("\\*\\(1e0/([a-zA-Z][a-zA-Z0-9_]*)\\)", "/\\1",
                      rate)  # Replace *1/x by /x to avoid dangling *1 in fractions
        rate = re.sub("\\+ *\\-", "-", rate)  # Replace + - by -

        rate = re.sub("\[", "{", rate)
        rate = re.sub("\]", "}", rate)

        # store for debugging
        rateTexAfterShortcutsReplaced = copy(rate)
        # transform to LaTeX format
        # keep trying
        while True:
            try:
                if opts.latex_backend == "pytexit":
                    rateTex = pytexit.for2tex(rate, print_latex=False, print_formula=False)
                    rateTex = rateTex[2:-2]
                elif opts.latex_backend == "sympy":
                    rateTex = sp.latex(eval(rate, symboltable))
                break
            # undefined variable will become a symbol
            except NameError as err:
                print("Name error in rate", err)
                varIssue = str(err).split("'")[1]
                symboltable[varIssue] = sp.Symbol(varIssue)

            # special case rate will be prited as it is
            except SyntaxError as err:
                if rate == "@xsecFile=SWRI": return "= \\textrm{table}", message, 1
                print("Syntax Error in rate", err)
                print("in ", rate)
                return "=" + rate, message

            # special case rate will be prited as it is
            except ValueError as err:
                print("Value Error in rate", err)
                return "=" + rate, message, 1

        # store for debugging
        rateTextAfterSympy = copy(rateTex)

        if opts.latex_backend == "pytexit":
            rateTex = utils.replaceSymbols(rateTex)

        # store for debugging
        rateTextAfterReplaceSymbols = copy(rateTex)

        # fix mistakes by sympy
        # no 10^{} for short rates
        # for t in Tsymbols:
        #    if t in rateTex:
        #        break
        #    else:
        #        print "YYYYYYYYYYYYYYYY"
        #        print rateTex
        #        rateTex = rateTex.replace("e-", "\\cdot 10^{-") + "}X"
        #        break

        # use a\cdot 10^-b for reaction that are just numbers
        try:
            float(rateTex.replace("=", ""))
            if "e-" in rateTex:
                rateTex = rateTex.replace("e-", "\\cdot 10^{-") + "}"
        except:
            pass

        # store for debugging
        rateTextAfterReplaceExpNot = copy(rateTex)

        # old algorithm
        #
        # #it automatically makes a fraction out of negative exponents
        # # TODO: make more generic, it fails with some reactions
        # stringFrac = r'\frac{1}{'
        # if stringFrac in rateTex: and not debug:
        #     for Tsym in Tsymbols:
        #         rateTex = rateTex.replace(stringFrac + Tsym + "^{", "{"+Tsym+"^{-")
        #     #change the order of the factors to match modified Arrhenius
        #     #avoid chaging stuff in more complex rates
        #     if len(rateTex) < maxRateLength:
        #         rateTexSplit = re.split("\}\}",rateTex)
        #         beta = rateTexSplit[0][1:]+"}" #beta containing factor
        #         if "exp" in rateTexSplit[-1]:
        #             alphaGamma = rateTexSplit[-1].split("\operatorname") #alpha and gamma factor
        #             rateTex = alphaGamma[0] + beta + "\operatorname" + alphaGamma[1]
        #         else:
        #             rateTex = rateTexSplit[-1] + beta

        # mistake: large fractions
        # solution: to the power -1 (solution can be improved)
        # TODO: make more generic, it fails with some reactions
        if rateTex.startswith(r"\frac") and not debug:
            cnt = 0
            pieces = []
            # get content between parenteses for each level
            parenticList = utils.getParentheticContents(rateTex, "{}")
            # only get the first \frac{}{} parts
            for part in parenticList:
                if part[0] == 0 and cnt < 3:
                    pieces.append(part[1])
                    cnt = cnt + 1

            # replace \frac{a}{b} with a*b^{-1}
            firstTerm = pieces[0]
            secondTerm = pieces[1]
            InvsecondTerm = None
            fracStringOriginal = r"\frac{" + firstTerm + "}{" + secondTerm + "}"
            restStringOriginal = rateTex.replace(fracStringOriginal, '')

            # invert denominator with changed sign
            if secondTerm.count('^{') == 1:
                InvsecondTerm = secondTerm.replace('^{', '^{-')
            # remove '1' as numerator
            if firstTerm.strip() == '1':
                firstTerm = ''

            if InvsecondTerm:
                fracStringReplace = firstTerm + InvsecondTerm
            # Put parenteses around first term if composite term
            elif " + " not in firstTerm or " - " not in firstTerm:
                fracStringReplace = firstTerm + "\left[" + secondTerm + r"\right]^{-1}"
            else:
                fracStringReplace = "\left[" + firstTerm + r"\right]\left[" + secondTerm + r"\right]^{-1}"

            # change the order of the factors to match modified Arrhenius
            # avoid changing stuff in more complex rates
            # TODO: Replace with regex expression similar to
            # utils.replaceFracsWithInvDenominators()
            if len(rateTex) < opts.max_fraction_length:
                # rateTexSplit = re.split("\}\}",rateTex)
                # print rateTexSplit
                # beta = rateTexSplit[0][1:]+"}" #beta containing factor

                if restStringOriginal.count('exp') == 1:
                    alphaFactor, gammaFactor = restStringOriginal.split("\operatorname{exp}")  # alpha and gamma factor
                    rateTex = alphaFactor + fracStringReplace + "\operatorname{exp}" + gammaFactor
                elif " + " in restStringOriginal:
                    splitted = restStringOriginal.split(" + ")
                    rateTex = splitted[0] + fracStringReplace + " + " + splitted[1]
                elif " - " in restStringOriginal:
                    splitted = restStringOriginal.split(" - ")
                    rateTex = splitted[0] + fracStringReplace + " - " + splitted[1]
                else:
                    rateTex = restStringOriginal + fracStringReplace

        # store for debugging
        rateTextAfterLargeFractsToExp = copy(rateTex)

        # replaces fractions with numerator = 1 and denominators with power
        # their inversed denominator
        # TODO: This functions can/should be used for more general purposes
        rateTex = utils.replaceFracsWithInvDenominators(rateTex, numerator='1')

        # remove unwanted zeros
        rateTex = re.sub(r"(\d+\.[1-9]*)0*(?=\D)", r"\1", rateTex)
        rateTex = re.sub(r"(\d+)\.(?=\D)", r"\1", rateTex)
        rateTex = re.sub(r"([^0-9\.])0*(\d+\.*)", r"\1\2", rateTex)

        # truncate numbers at 5 decimal places
        # TODO: properly round the values
        rateTex = re.sub(r'(\d+\.[0-9]{' + str(opts.truncate_numbers) + '})\d*', r'\1', rateTex)

        rateTex = re.sub("([ \)_]*)idx_{([A-Za-z0-9_]{1,})}", r"\1idx_\2", rateTex)
        rateTexIdxReplaced = copy(rateTex)

        # Rename number density, e.g. n(idx_H) -> n_idx_H
        # "(?<!\l)" is a negative lookback that ensures "\ln{...}" is not
        # captured by this regex
        rateTex = re.sub(r"([^A-Za-z])?(?<!\l)n\{(\\left|) *\( *([A-Za-z0-9_]{1,}) *(\\right|) *\)\}", r"\1n_{\3}",
                         rateTex)
        # Rename species id to name wrapped in \ch, e.g. idx_H2 -> \ch{H2}
        rateTex = re.sub("([ \)_]*)idx_([A-Za-z0-9_]{1,})( *)", r"\1\ch{\2}\3", rateTex)

        rateTexBeforeReplacePm = copy(rateTex)
        rateTex = re.sub("\\+ *\\-", "-", rateTex)  # Replace + - by -

        # store rate before breaking
        rateTextFull = copy(rateTex)

        # break long rates in multiple lines
        if len(rateTex) > opts.max_fraction_length:
            rateTex, numlines = self.breakRateTex(rateTex)
        else:
            # add LaTeX symbols
            rateTex = " = " + rateTex
            numlines = 1

        message += "\n%These comments below are for debugging, ignore them"
        message += "\n%original rate: " + originalRate
        message += "\n%sympy friendly rate: " + rate
        message += "\n%after shortcuts replacing: " + rateTexAfterShortcutsReplaced
        message += "\n%rate after sympy: " + rateTextAfterSympy
        message += "\n%rate after replace symbols" + rateTextAfterReplaceSymbols
        message += "\n%rate after replace exp. not." + rateTextAfterReplaceExpNot
        message += "\n%rate after replace large fracts. w. exp. not" + rateTextAfterLargeFractsToExp
        message += "\n%rate after replacing idx: " + rateTexIdxReplaced
        message += "\n%rate before replaceing +-: " + rateTexBeforeReplacePm
        message += "\n%full rate (before breaking): " + rateTextFull

        return rateTex, message, numlines

    # ****************
    # look up auto rate in database
    def find_auto_rate(self, reactants, products, reactionType):
        import os, re

        reactants = ",".join(map(lambda x: x.name.replace("+", '\+'), reactants))
        products = " *, *".join(map(lambda x: x.name.replace("+", '\+'), products))

        expr_w_limits = """ *@type: *{}ion
 *@reacts: *{}
 *@prods: *{}
 *@limits: ([0-9\.de\+\<\>]*), *([0-9\.de\+\<\>]*)
 *@rate: *(.*)""".format(reactionType, reactants, products)

        expr_wo_limits = """ *@type: *{}ion
 *@reacts: *{}
 *@prods: *{}
 *@rate: *(.*)""".format(reactionType, reactants, products)

        re_w_limits = re.compile(expr_w_limits, re.MULTILINE)
        re_wo_limits = re.compile(expr_wo_limits, re.MULTILINE)

        self.database_path = "../../data/database"
        for filename in os.listdir(self.database_path):
            fullpath = os.path.join(self.database_path, filename)
            if os.path.isfile(fullpath):
                with open(fullpath, 'r') as f:
                    contents = f.read()
                    result = re_wo_limits.search(contents) or re_w_limits.search(contents)
                    if (result):
                        groups = list(result.groups())
                        rate = groups.pop()
                        if len(groups) > 0:
                            Emax = groups.pop()
                        if len(groups) > 0:
                            Emin = groups.pop()
                        return (rate, [Emin.replace("d", "e")], [Emax.replace("d", "e")])

        # raise Exception("Auto rate not found in database for reaction " + self.getVerbatim())
        print("Auto rate not found in database for reaction " + self.getVerbatim())
        return 'auto', None, None

    # ****************
    # break long rates and put in LateX table format
    def breakRateTex(self, rate):
        from options import latexoptions as opts

        rate = utils.raiseBracketsOnOperators(rate)
        rate = utils.replaceLongFracByDivide(rate, maxlen=opts.max_fraction_length)
        rate = utils.replaceLongExponentByHat(rate, maxlen=opts.max_exponent_length)
        rate = utils.replaceLeftRightbyBigLR(rate)
        rate, numlines = utils.breakLatexEquation(rate, maxlen=opts.max_fraction_length)
        rate = "\\begin{aligned}[t] = &" + rate + "\\end{aligned}"
        return rate, numlines

    # ****************
    # break long rates and put in LateX table format
    def breakRateTexDeprecated(self, rate):
        nextReplace = False

        # get rid of unclosed "{}" on a line, when breaking equation
        parList = utils.getParentheticContents(rate, "{}")
        for part in parList:
            # replace when level 0 and after exp
            if part[0] == 0 and nextReplace:
                rate = rate.replace("{" + part[1] + "}", part[1])
                nextReplace = False
            # only for exp
            if part == (0, "exp"):
                # next level 0 parenteses need to be replaced
                nextReplace = True

        rate = rate[1:-1]  # remove $ signs
        rate = rate.replace(" + ", " \\\\ \n& + ")
        rate = rate.replace(" - ", " \\\\ \n& - ")
        rate = "\\begin{aligned}[t] & " + rate + "\\end{aligned}$"
        warning = ""

        newLines = []
        for line_rn in rate.split("&"):
            ends_rn = line_rn.strip().endswith("\\\\")
            line = line_rn.strip().replace("\\\\", "")
            Nleft = line.count("\\left")
            Nright = line.count("\\right")
            if Nleft > Nright:
                warning = "\n%%*********************\n%% added some \\right to balance \\left"
                # auto replace was not always succesful...
                # newline = line.replace("\\left","\\left\\right.", Nleft-Nright)
                line += "\\right." * (Nleft - Nright)
                # rate = rate.replace(line, newline)
            elif Nleft < Nright:
                warning = "\n%%*********************\n%% added some \\left to balance \\right"
                # newline = line.replace("\\right","\\left.\\right", Nright-Nleft)
                line = "\\left." * (Nright - Nleft) + line
                # rate = rate.replace(line, newline)
            if ends_rn:
                line += "\\\\\n"
            newLines.append(line)

        rate = " & ".join(newLines)

        return rate, warning

    # ****************
    # temperature rage to LateX format
    def tempRange2latex(self, idxMerged):
        import sympy as sp
        num = sp.__version__.count('.') - 1
        sp_version = float(sp.__version__.rsplit('.', num)[0])
        if sp_version >= 1.3:
            print("ERROR: The LaTeX conversion currently only works with and older"
                  " version of SymPy (<1.3). Symbols no longer automatically"
                  " convert to functions when called."
                  )
            sys.exit()

        # change limits to uniform format
        low = utils.simplifyLimits(self.Tmin[idxMerged])
        high = utils.simplifyLimits(self.Tmax[idxMerged])

        # algorithm to account for differnt KROME formater of limits
        # e.g. with or without ">", "<",
        # put in correct order and switch symbols if needed
        if low:
            if ">" not in low:
                low = " > " + low
            else:
                low = low.replace(">", " > ").replace("> =", " >= ")
            if not (high):
                lowhigh = " T " + low + " K "
            else:
                low = low.split()[-1] + " " + "".join(low.split()[:-1])
                low = low.replace(">", "<")

        if high:
            if "<" not in high:
                high = " < " + high
            else:
                high = high.replace("<", " < ").replace("< =", " <= ")
            if not low:
                low = ""
            lowhigh = low + " T " + high + " K "

        if not low and not high:
            lowhigh = ""

        # change limits symbols to latex format
        lowhighTex = utils.limits2latex(lowhigh)

        # turn numbers into integers in latex format
        lowhighTex = [str(utils.char2int(part)) for part in lowhighTex.split()]
        limitTex = " ".join(lowhighTex)
        limitTex = utils.exp2latex(limitTex)

        return limitTex
