import sys, glob, os, shutil, hashlib
import reaction, utils, options, copy


class network:

    # ************************
    # class constructor
    def __init__(self, myOptions):

        self.species = []
        self.speciesDictionary = []
        self.myOptions = myOptions

        basePath = os.path.join(os.path.dirname(__file__), "..")
        basePath = os.path.abspath(basePath)

        # get network file name
        fileName = myOptions.network

        # read atoms
        absPath = os.path.join(basePath, "atomlist.dat")
        atomSet = utils.getAtomSet(absPath)
        self.atomSet = atomSet

        # load thermochemical data
        absPath = os.path.join(basePath, "thermo30.dat")
        self.thermochemicalData = utils.getThermochemicalData(absPath)

        # load polarizability data
        absPath = os.path.join(basePath, "polarizability.dat")
        self.polarizabilityData = utils.getPolarizabilityData(absPath)

        # read shortcuts
        shortcuts = utils.getShortcuts()

        # read network file in KROME format
        if self.detectNetworkType(fileName) == "KIDA":
            self.readFileKIDA(fileName, atomSet)
        elif self.detectNetworkType(fileName) == "UMIST":
            self.readFileUMIST(fileName, atomSet)
        else:
            self.readFileKROME(fileName, atomSet, shortcuts)

        # merge reactions with multiple rates
        print("merging multiple reactions")
        self.mergeReactions()

        # get a list with missing reactions (reactants)
        self.missingReactions = self.getMissing()

        # get ranges from option file
        varRanges = dict()
        for rng in self.myOptions.range:
            # store range name and range limits
            (rangeName, rangeValue) = [x.strip() for x in rng.split("=")]
            varRanges[rangeName] = [float(x) for x in rangeValue.split(",")]
            print(rangeName, varRanges[rangeName])

        # clear folders (html and png)
        self.clearFolders()
        self.backupEvaluationJSON()

        # prepare subnetwork
        self.subNetwork(self.myOptions)

        # check missing xsecs file
        self.reportXsecsFiles()

        # loop on reactions to update attributes
        for myReaction in self.reactions:
            myReaction.updateAttributes()

        # loop on reactions to evaluate
        for myReaction in self.reactions:
            myReaction.evaluateRate(shortcuts, varRanges)

    # *****************************
    # assign index to species from info.log file produced by KROME
    def assignIndex(self, fileName="info.log"):

        inBlock = False  # flag in reading block
        self.indexDictionary = dict()  # index dictionary, key=species name, value=index

        # loop on log file
        for row in open(fileName):
            # check beginning line to read
            if row.startswith("# Species list with their indexes"):
                inBlock = True
                continue
            # if outside reading block skip
            if not inBlock:
                continue
            # check ending of reading block and break from reading
            if inBlock and row.strip() == "":
                break
            # replace tabs with spaces
            srow = row.strip().replace("\t", " ")
            # read data
            (idx, speciesName, idxF90) = [x.strip() for x in srow.split(" ") if x != ""]
            # assign data to dictionary
            self.indexDictionary[speciesName] = int(idx)

        # assign index to class species dictionary
        for species in self.species:
            # check for capital names
            if not species.name in self.indexDictionary:
                species.idx = self.indexDictionary[species.name.upper()]
            else:
                species.idx = self.indexDictionary[species.name]

    # ******************************
    # prepare complete documentation
    def makeAllDoc(self):

        # prepare graphs
        self.makeGraph()

        # prepare html
        self.makeHTML()

        # plot all rates
        self.plotRates()

        # make latex table of network
        # self.network2latex()

    # ********************
    # prepare all the html pages
    def makeHTML(self):

        # prepare html pages for species
        for mySpecies in self.getSpecies():
            mySpecies.makeHtmlPage(self)
            mySpecies.makeAllRatesHtmlPage(self, self.myOptions)

        # loop on reactions to evaluate
        for myReaction in self.reactions:
            myReaction.makeHtmlPage(self.myOptions, self)

        # create HTML pages
        self.makeHtmlIndex()
        self.makeHtmlMissingBranchesIndex(self.myOptions)
        self.makeHtmlMissingReactionIndex(self.myOptions)
        self.makeHtmlReactionIndex(self.myOptions)
        self.makeHtmlSpeciesIndex()
        self.makeHtmlGraphIndex()
        self.makeHtmlAllRates(self.myOptions)
        self.makeHtmlMultipleRates(self.myOptions)

    # **********************
    # plot all rates to PNG
    def plotRates(self):

        # delete plots that are not changed
        self.deleteChangedPNGs(self.myOptions)

        # plotting rates
        print("plotting rates...")
        icount = 0
        # loop on reactions to plot
        for myReaction in self.reactions:
            print(str(int(icount * 1e2 / len(self.reactions))) + "%", myReaction.getVerbatim())
            myReaction.plotRate(self.myOptions)
            icount += 1

    # ********************
    # make a LaTeX table of the network
    def network2latex(self, networkLatex="NetworkLatex.tex"):
        from options import latexoptions as opts
        # NOTE: Make sure the network input file has incrementing reaction indices.
        # If double indices exist, the LaTeX table will be incorrect.

        # list with all temperature shortcuts element = (var, replaceWith)
        shortcutsTemperature = utils.getShortcutsLatex()
        deferredShortcuts = utils.getDeferredShortcuts()
        cntMergedReactions = 0
        cntTotalReactions = 1
        cntAllReactions = 0
        linesOnPage = 0
        replace_references = {}
        ref_cnt = 1

        if opts.sorted_by == 'alphabetic':
            def reaction_sort_function(x):
                return [xx.name for xx in x.reactants]
        else:
            # default is on index
            def reaction_sort_function(x):
                return x.uid

        with open(networkLatex, "w") as fileOutput:
            # dump header of the file
            self.dumpLatexTableHeader(fileOutput)
            # dump beginning of table float, table environment and header row
            self.dumpLatexBeginTable(fileOutput)
            # loop on reactions to evaluate
            for myReaction in sorted(self.reactions, key=reaction_sort_function):
                # list with all variable shortcuts (excl. temperature ones)
                shortcutsVariables = myReaction.shortcuts

                nrates = len(myReaction.rate)
                rateColumns = []
                messages = []
                linesInReaction = 0
                for cnt in range(nrates):
                    latexColumns, message, linesInRate = myReaction.reaction2latex(shortcutsTemperature,
                                                                                   shortcutsVariables,
                                                                                   deferredShortcuts,
                                                                                   cntMergedReactions,
                                                                                   cnt,
                                                                                   cntTotalReactions,
                                                                                   cntAllReactions)

                    # Make sure the references are in ascending order
                    # This needs to be done if the reactions are not sorted
                    # on index
                    ref = latexColumns[-1]
                    if ref == '\\dbDefault{}':
                        pass
                    else:
                        try:
                            new_ref = replace_references[ref]
                        except KeyError:
                            replace_references[ref] = str(ref_cnt)
                            ref_cnt += 1
                            new_ref = str(ref_cnt - 1)
                        latexColumns[-1] = new_ref

                    rateColumns.append(latexColumns)
                    messages.append(message)
                    linesInReaction += linesInRate
                linesOnPage += linesInReaction

                if linesOnPage > opts.lines_per_page:
                    # Dump page break
                    self.dumpLatexEndTable(fileOutput)
                    self.dumpLatexBeginTable(fileOutput, first=False)
                    linesOnPage = linesInReaction

                for columns, message in zip(rateColumns, messages):
                    if message:
                        fileOutput.write(message + "\n")
                    self.dumpLatexTable(columns, fileOutput)

                cntTotalReactions += 1
            # Dump end of table environment and - float
            self.dumpLatexEndTable(fileOutput, last=True)

        self.dumpDeferredShortcuts(shortcutsTemperature, shortcutsVariables, deferredShortcuts)
        # update referenceId with new ascending reference indices
        self.referenceId = {key: int(replace_references[str(val)]) for key, val in self.referenceId.items()}
        self.dumpLatexReferences()

    # ****************
    # dump colums in LateX table format
    def dumpLatexTableHeader(self, tableFile):

        tableFile.write("%********************************\n")
        tableFile.write("% Include \\usepackage{chemformula}\n")
        tableFile.write("% Chemical reaction network (LaTeX format)\n")
        tableFile.write("% Table columns format {" + 5 * "l" + "}\n")
        tableFile.write("% Set \"h\" as table-spec for the column to hide it\n")
        tableFile.write("%%\\newcolumntype{h}{>{\setbox0=\hbox\\bgroup}c<{\egroup}@{}}\n")

    # ****************
    # dump headers for table float and - environment
    def dumpLatexBeginTable(self, tableFile, first=True):
        arg = "\\chemicalNetworkTableCaption"
        if not first: arg += "Cont"
        if first: arg += "\n\\label{\\chemicalNetworkTableLabel}"
        tableFile.write("\\begin{chemical_network_table}{" + arg + "}\n")
        tableFile.write("  \\begin{chemical_network_tabular}\n")

    # ****************
    # dump footers for table float and - environment
    def dumpLatexEndTable(self, tableFile, last=False):
        tableFile.write("  \\end{chemical_network_tabular}\n")
        if last:
            tableFile.write("  \\chemicalNetworkTableNotes\n")
        tableFile.write("\\end{chemical_network_table}\n\n")

    # ****************
    # build colums in LateX table format
    def dumpLatexTable(self, columns, tableFile):
        row = " & ".join(columns) + " \\\\"
        tableFile.write(row + "\n")

    # ****************
    # dump a list of equations for deferred variables
    def dumpDeferredShortcuts(self, temperatureShortcuts, variableShortcuts, deferredShortcuts,
                              filename="NetworkLatexSymbols.tex"):
        import re
        from options import latexoptions as opt
        import utils
        if opt.latex_backend == "pytexit":
            import pytexit
        elif opt.latex_backend == "sympy":
            import sympy as sp
            num = sp.__version__.count('.') - 1
            sp_version = float(sp.__version__.rsplit('.', num)[0])
            if sp_version >= 1.3:
                print("ERROR: The LaTeX conversion currently only works with and older"
                      " version of SymPy (<1.3). Symbols no longer automatically"
                      " convert to functions when called."
                      )
                sys.exit()

        symboltable = utils.getSymbolTable()
        with open(filename, "w") as fileOutput:
            for expr, symbol in deferredShortcuts.iteritems():
                expr = utils.replaceShortcuts(expr, variableShortcuts, [])
                expr = utils.replaceShortcuts(expr, temperatureShortcuts[1:], [])
                expr = expr.replace("dexp", "exp")
                expr = expr.replace("d", "e")
                expr = expr.replace("log", "ln")
                expr = expr.replace("ln10", "log")
                expr = expr.replace("_8", "")
                expr = re.sub("\\*1e0/([a-zA-Z][a-zA-Z0-9_]*)", "/\\1",
                              expr)  # Replace *1/x by /x to avoid dangling *1 in fractions
                expr = re.sub("\\+ *\\-", "-", expr)  # Replace + - by -
                if opt.latex_backend == "pytexit":
                    exprTex = pytexit.for2tex(expr, print_latex=False, print_formula=False)
                    exprTex = exprTex[2:-2]
                else:
                    while True:
                        try:
                            exprTex = sp.latex(eval(expr, symboltable))
                            break
                        # undefined variable will become a symbol
                        except NameError as err:
                            print("Name error in expression", err)
                            varIssue = str(err).split("'")[1]
                            symboltable[varIssue] = varIssue

                # use a\cdot 10^-b for reaction that are just numbers
                try:
                    float(exprTex.replace("=", ""))
                    exprTex = exprTex.replace("e-", "\\cdot 10^{-") + "}"
                except:
                    pass

                # remove unwanted zeros
                exprTex = re.sub(r"(\d+\.[1-9]*)0*(?=\D)", r"\1", exprTex)
                exprTex = re.sub(r"(\d+)\.(?=\D)", r"\1", exprTex)
                exprTex = re.sub(r"0*(\d+\.*)", r"\1", exprTex)

                # pytexit doesn't replace symbols, so it is done here
                if opt.latex_backend == "pytexit":
                    exprTex = utils.replaceSymbols(exprTex)

                fileOutput.write("$" + symbol + "$ & $ = " + exprTex + "$ \\\\ \n")

    # ****************
    # make a Latex list of references for the reaction table
    def dumpLatexReferences(self, filename="NetworkLatexReferences.tex"):
        with open(filename, "w") as fileOutput:
            refs = [(id, ref) for ref, id in self.referenceId.items()]
            refs = ["(" + str(id) + ") " + ref for id, ref in sorted(refs)]
            if len(refs) > 1:
                refstex = ", ".join(refs[:-1]) + " and " + refs[-1]
            elif len(refs) == 1:
                refstex = refs[0]
            else:
                refstex = ""
            fileOutput.write(refstex)

    # ********************
    # check xsec folder and write xsec status to log file
    def reportXsecsFiles(self, logName="xsecs.log", folder="xsecs/"):

        fout = open(logName, "w")
        # loop on species
        for spec in self.species:
            fname = folder + spec.name + ".dat"
            status = "missing <<<"
            # check if file is present
            if os.path.exists(fname):
                status = "present"
            specName = spec.name + (" " * (30 - len(spec.name)))
            fout.write(specName + status + "\n")
        fout.close()
        print("xsecs data report written to " + logName)

    # ********************
    def detectNetworkType(self, fileName):
        fh = open(fileName)
        for row in fh:
            srow = row.strip()
            if srow == "": continue
            if srow.startswith("#"): continue
            if "," in srow: return "KROME"
            if srow.count(":") > 5: return "UMIST"
        fh.close()
        return "KIDA"

    # ************************
    def readFileKROME(self, fileName, atomSet, shortcuts):

        # default format
        reactionFormat = "@format:idx,R,R,R,P,P,P,P,Tmin,Tmax,rate"
        shortcutsList = []  # used in network2latex, list because dict loses order
        inBlockCR = False
        inBlockPhoto = False
        inBlockCatalysis = False

        self.reactions = []
        reactionId = {}
        self.referenceId = {}
        nextReactionId = 0
        nextReferenceId = 0
        print("reading network " + fileName)
        # open file to read network
        fh = open(fileName)
        reference = None
        for row in fh:
            srow = row.strip()
            if srow == "": continue
            if srow.startswith("#@ref:"):
                reference = srow.replace("#@ref:", "")
                if not reference in self.referenceId:
                    nextReferenceId += 1
                    self.referenceId[reference] = nextReferenceId
            if srow.startswith("#"): continue

            # get format if token found
            if srow.startswith("@format:"):
                reactionFormat = srow
                continue

            # get extra variable definition if token found
            if srow.startswith("@var:"):
                (variable, expression) = [x.strip() for x in srow.replace("@var:", "").split("=")]
                shortcuts[variable] = expression
                # check if not already in temperature shortcuts
                if not (utils.isTemperatureShortcut(variable)):
                    shortcutsList.append((variable, expression))

            if srow.lower().startswith("@cr_start") or srow.lower().startswith("@cr_begin"):
                inBlockCR = True
            if srow.lower().startswith("@cr_stop") or srow.lower().startswith("@cr_end"):
                inBlockCR = False

            if srow.lower().startswith("@photo_start") or srow.lower().startswith("@photo_begin"):
                inBlockPhoto = True
            if srow.lower().startswith("@photo_stop") or srow.lower().startswith("@photo_end"):
                inBlockPhoto = False

            if srow.lower().startswith("@catalysis_start") or srow.lower().startswith("@catalysis_begin"):
                inBlockCatalysis = True
            if srow.lower().startswith("@catalysis_stop") or srow.lower().startswith("@catalysis_end"):
                inBlockCatalysis = False

            # change reaction type to CR
            reactionType = "standard"  # default
            if inBlockCR: reactionType = "CR"
            if inBlockPhoto: reactionType = "photo"
            if inBlockCatalysis: reactionType = "catalysis"

            # skip other tokens
            if srow.startswith("@"): continue

            # parse row line for reaction
            if reference is None:
                ref = None
            else:
                ref = (self.referenceId[reference], reference)

            myReaction = reaction.reaction(srow, reactionFormat, atomSet, reactionType, self.species,
                                           shortcutsList, reference=ref)

            # update reaction id
            key = (str(myReaction.reactants), str(myReaction.products), str(myReaction.reactionType))
            if key in reactionId:
                myReaction.uid = reactionId[key]
            else:
                nextReactionId += 1
                myReaction.uid = nextReactionId
                reactionId[key] = myReaction.uid

            # add parsed reaction to reactions structure in network
            self.reactions.append(myReaction)

        fh.close()

    # **************
    def readFileKIDA(self, fileName, atomSet):
        self.reactions = []
        print("reading network " + fileName)

        fh = open(fileName)
        for row in fh:
            srow = row.strip()
            if srow == "": continue
            if srow.startswith("#"): continue
            # skip kida comment
            if srow.startswith("!"): continue

            reactionType = "KIDA"
            reactionFormat = ""

            # parse row line for reaction
            myReaction = reaction.reaction(srow, reactionFormat, atomSet, reactionType, self.species)

            # add parsed reaction to reactions structure in network
            self.reactions.append(myReaction)

        fh.close()

    # **************
    def readFileUMIST(self, fileName, atomSet):
        self.reactions = []
        print("reading network " + fileName)

        fh = open(fileName)
        for row in fh:
            srow = row.strip()
            if srow == "": continue
            if srow.startswith("#"): continue

            reactionType = "UMIST"
            reactionFormat = ""

            # parse row line for reaction
            myReaction = reaction.reaction(srow, reactionFormat, atomSet, reactionType, self.species)

            # add parsed reaction to reactions structure in network
            self.reactions.append(myReaction)

        fh.close()

    # **************
    # get all network species
    def getSpecies(self):
        return self.species

    # **************
    # get all network species in a dictionary with NAME as keys
    def getSpeciesDictionary(self):
        if self.speciesDictionary != []:
            return self.speciesDictionary

        # unique list of species
        self.speciesDictionary = dict()
        for mySpecies in self.getSpecies():
            self.speciesDictionary[mySpecies.name] = mySpecies

        return self.speciesDictionary

    # *****************
    # get species object from species list
    def getSpeciesByName(self, name):
        speciesDictionary = self.getSpeciesDictionary()
        if name not in speciesDictionary:
            print("ERROR: species " + name + " unknown!")
            print(" Available species are:")
            print(", ".join(speciesDictionary.keys()))
            sys.exit()
        return speciesDictionary[name]

    # *******************
    def getAtoms(self):
        atoms = []
        for mySpecies in self.getSpecies():
            atoms += mySpecies.atoms
        return list(set(atoms))

    # **************
    # clean temporary folders
    def clearFolders(self):

        # files to be deleted (folder:extension)
        # png and json are cleaned later if different
        folders = {"htmls": "html",
                   "epss": "eps",
                   "dots": "dot"}
        # "pngs":"png", \
        # "evals":"json", \

        # loop on folders and extensions
        for (path, extension) in folders.items():
            removePath = path + "/*." + extension
            print("removing " + removePath)
            # get files list
            filelist = glob.glob(removePath)
            # remove files
            for fname in filelist:
                os.remove(fname)

    # **************
    # backup json files with rate evaluations (to avoid replotting)
    def backupEvaluationJSON(self):

        backupPath = "evals/*.json"
        # get files list
        filelist = glob.glob(backupPath)
        # copy files
        for fname in filelist:
            shutil.copyfile(fname, fname + ".bak")

    # **************
    # delete PNG files with no json.bak file or different MD5
    # after this remove all json.bak files
    def deleteChangedPNGs(self, myOptions):

        # get file lists
        pngList = glob.glob("pngs/rate_*.png")

        # loop on variable ranges
        for rng in myOptions.range:
            # get range name
            rngName = rng.split("=")[0].strip()
            # produce end of the file name
            fileNameEnd = "_" + rngName + ".png"
            # copy files
            for pngFname in pngList:
                survived = True
                # if png files end with range name
                if fileNameEnd in pngFname:
                    # create json and json.bak filenames
                    jsonFname = pngFname.replace(fileNameEnd, "").replace("pngs/", "evals/") + ".json"
                    jsonBakFname = jsonFname + ".bak"
                    # if bak file is missing remove png file
                    if not os.path.exists(jsonBakFname):
                        # print "missing json bak", pngFname
                        survived = False
                        os.remove(pngFname)
                    else:
                        # check MD5 of json and json.bak
                        md5 = hashlib.md5(open(jsonFname).read()).hexdigest()
                        md5Bak = hashlib.md5(open(jsonBakFname).read()).hexdigest()
                        # if md5 are different remove png file
                        if md5 != md5Bak:
                            # print "diff MD5", pngFname,md5, md5Bak
                            survived = False
                            os.remove(pngFname)
                    # if(survived): print "survived: "+pngFname

        # remove json.bak files
        bakList = glob.glob("evals/*.bak")
        for jsonBak in bakList:
            os.remove(jsonBak)

    # **************
    # merge reactions with multiple rates
    def mergeReactions(self):
        reactionList = dict()
        # merge reactions
        for myReaction in self.reactions:
            myHash = myReaction.getReactionHash()
            if myHash in reactionList:
                reactionList[myHash].merge(myReaction)
            else:
                reactionList[myHash] = myReaction
        unsortedReactions = [v for (k, v) in reactionList.items()]
        self.reactions = sorted(unsortedReactions, key=lambda x: x.index)

    # ***************
    def makeGraph(self):

        import subprocess

        # maximum number of nodes
        maxNodes = 50

        # loop on all available atoms in the network
        for refAtom in self.getAtoms():
            print("creating graphs for " + refAtom + "...")

            # connection dictionary
            networkMap = dict()
            # loop on reactions
            for myReaction in self.reactions:
                # loop on reactants
                for myR in myReaction.reactants:
                    # skip species if reference atom is not contained
                    if (refAtom not in myR.atoms) and (refAtom is not None): continue
                    # get reactant partners in this reaction (other reactants)
                    partners = [x.name for x in myReaction.reactants if (x.name != myR.name)]
                    # loop on products
                    for myP in myReaction.products:
                        if myP.name == myR.name: continue
                        # skip species if reference atom is not contained
                        if (refAtom not in myP.atoms) and (refAtom is not None): continue
                        # create dictionary if not present
                        if myR.name not in networkMap:
                            networkMap[myR.name] = dict()
                        # create list if not present
                        if myP.name not in networkMap[myR.name]:
                            networkMap[myR.name][myP.name] = []
                        # add partners to the matrix
                        networkMap[myR.name][myP.name] += partners

            # dot, png, eps files
            fname = "dots/network_atom_" + str(refAtom) + ".dot"
            fnamePNG = "pngs/network_atom_" + str(refAtom) + ".png"
            fnameEPS = "epss/network_atom_" + str(refAtom) + ".eps"

            # write dot
            fout = open(fname, "w")
            fout.write("digraph G{\n")

            if len(networkMap.keys()) == 0: continue

            # skip network with too many nodes
            if len(networkMap.keys()) > maxNodes:
                print("WARNING: too many nodes, skipping graph!")
                continue

            nodeLabels = dict()
            # loop on starting edges
            for myR in networkMap.keys():
                nodeLabels[myR.lower()] = myR
                # loop on ending edges (and partners)
                for (myP, partners) in networkMap[myR].items():
                    nodeLabels[myP.lower()] = myP
                    labelList = sorted(list(set(partners)))
                    label = ",".join(labelList)
                    # write edge to dot
                    fout.write("\"" + myR.lower() + "\" -> \"" + myP.lower() + "\" [label=\"" \
                               + label + "\"];\n")

            for (k, v) in nodeLabels.items():
                fout.write("\"" + k + "\" [label=\"" + v + "\"];\n")

            fout.write("}\n")
            fout.close()

            # write PNG with dot
            myCall = "dot " + fname + " -Tpng -o " + fnamePNG
            subprocess.call(myCall.split(" "))

            # write EPS with dot
            myCall = "dot " + fname + " -Teps -o " + fnameEPS
            subprocess.call(myCall.split(" "))

    # ****************
    # return all the possible reactants combinations with the current species
    # as an array, including unimol and bimol reactions
    def getAllReactantsCombinations(self, returnNames=True, includeRepulsive=False):

        # get list of species
        species = self.getSpecies()

        # unimolecular reactants to be skipped
        skipUnimol = ["CR", "E"]

        allReactants = []
        # loop on species (reactant 1)
        for species1 in species:
            # add unimolecular (only neutral or anion)
            if (species1.charge < 1) and (species1.name not in skipUnimol):
                if returnNames:
                    allReactants.append([species1.name])
                else:
                    allReactants.append([species1])
            # loop on species (reactant 2) for bimolecular
            for species2 in species:
                # exclude repulsive
                if (species1.charge * species2.charge > 0) and not includeRepulsive: continue

                # exclude double ionization CR
                hasCR = (species1.name == "CR" or species2.name == "CR")
                hasElectron = (species1.name == "E" or species2.name == "E")
                hasCation = (species1.charge > 0 or species2.charge > 0)
                if hasCation and hasCR: continue
                if hasCR and hasElectron: continue

                # add bimolecular
                if returnNames:
                    reactants = [species1.name, species2.name]
                    allReactants.append(sorted(reactants))
                else:
                    allReactants.append([species1, species2])

        return allReactants

    # ****************
    # get a list of the missing reactions (reactants only)
    def getMissing(self):

        # get list of species
        species = self.getSpecies()

        # limit to produce missing reactions
        expectedReactionsMax = int(1e4)

        # evaluate expected reactions
        expectedReactions = len(species) ** 2 + len(species)
        print("expected reactions:", expectedReactions)

        # skip if too many reactions
        if expectedReactions > expectedReactionsMax:
            print("WARNING: too many expected reactions skipping missing reactions!")
            return []

        # get all possible reactants combinations
        allReactants = self.getAllReactantsCombinations()

        # get species dictionary to return species objects instead of names
        speciesDictionary = self.getSpeciesDictionary()

        # create a dictionary of products using exploded as key
        # e.g. -_C_H_H_O for H2O + C-
        productDictionary = dict()
        # loop on all reactants (products are the same since species1+species2 generated)
        for reactants in allReactants:
            reactExploded = []
            # loop on reactants for exploded (not explodedFull contains charge signs)
            for RR in reactants:
                reactExploded += speciesDictionary[RR].explodedFull
            # concatenate exploded, e.g. -_C_H_H_O for H2O + C-
            reactExp = ("_".join(sorted(reactExploded)))
            # remove mutual sign neutralization
            reactExp = reactExp.replace("+_-_", "")
            # store products in a list (init list if not present)
            if reactExp not in productDictionary:
                productDictionary[reactExp] = []
            # store unique products
            if reactants not in productDictionary[reactExp]:
                productDictionary[reactExp].append(reactants)

        # write number of potential reactions found
        print("expected reactions (after criteria):", len(allReactants))

        # store reactants for each reaction
        networkReactants = []
        for myReaction in self.reactions:
            reactants = sorted([x.name for x in myReaction.reactants])
            networkReactants.append(reactants)

        # search for missing rates
        missing = []
        # loop on all possible reactants groups
        for reactant in allReactants:
            found = False
            # loop on network reactant groups
            for rectantNtw in networkReactants:
                if reactant == rectantNtw:
                    found = True
                    break
            # store not found if not present already
            if not found and reactant not in missing: missing.append(reactant)

        # sort missing reactions by reactants name
        missing = sorted(missing, key=lambda x: ("_".join(x)))
        print("missing reactions:", len(missing))

        missingReactions = []
        # loop on missing reactions to search for possible branches
        for reactants in missing:
            # create exploded reaction to compare and get branches
            reactExploded = []
            for RR in reactants:
                reactExploded += speciesDictionary[RR].explodedFull
            reactExp = ("_".join(sorted(reactExploded)))
            # remove mutual sign neutralization
            reactExp = reactExp.replace("+_-_", "")
            # skip reactants==products, e.g. H+OH -> H+OH (and convert names into species objects)
            productsList = [[speciesDictionary[x] for x in products] \
                            for products in productDictionary[reactExp] if (products != reactants)]
            productsListPrune = []
            for products in productsList:
                charges = [x.charge for x in products]
                # skip reaction that generates anion-cation products
                if (max(charges) + min(charges) == 0) and max(charges) > 0: continue
                productsListPrune.append(products)

            # convert names in to species objects
            reatantsList = [speciesDictionary[x] for x in reactants]
            missingReactions.append({"reactants": reatantsList, "products": productsListPrune})

        # return list of objects reactants
        return missingReactions

    # ****************
    # given a list of species objects return the sum of their exploded
    # taking care of charges
    def getExplodedSpecies(self, speciesList):

        rexp = []
        # sum exploded
        for species in speciesList:
            rexp += species.explodedFull
        # wrap into underscores
        rexpFull = "_" + ("_".join(sorted(rexp))) + "_"
        # remove +/- one by one until both are present
        while ("+" in rexpFull) and ("-" in rexpFull):
            rexpFull = rexpFull.replace("+", "", 1)
            rexpFull = rexpFull.replace("-", "", 1)
        # remove +/E one by one until both are present
        while ("_+_" in rexpFull) and ("_E_" in rexpFull):
            rexpFull = rexpFull.replace("+", "", 1)
            rexpFull = rexpFull.replace("E", "", 1)

        # for this hash E is just represented  as -
        rexpFull = rexpFull.replace("_E_", "_-_")
        # remove CRs
        rexpFull = rexpFull.replace("CR", "")
        # remove double underscores
        rexpFull = rexpFull.replace("__", "_")

        # sort again replaced
        rexpFull = ("_".join(sorted(rexpFull.split("_"))))

        # remove leading and trailing underscores
        while rexpFull.startswith("_"):
            rexpFull = rexpFull[1:]
        while rexpFull.endswith("_"):
            rexpFull = rexpFull[:-1]

        return rexpFull

    # ****************
    def getMissingBranch(self):

        missingBranch = []

        species = self.getSpecies()

        # get all possible species combinations
        allSpecies = dict()
        for spec1 in species:
            unimol = [spec1]
            rHash = self.getExplodedSpecies(unimol)
            if rHash not in allSpecies: allSpecies[rHash] = []
            allSpecies[rHash].append(unimol)
            # bimolecular
            for spec2 in species:
                bimol = [spec1, spec2]
                rHash = self.getExplodedSpecies(bimol)
                if rHash not in allSpecies: allSpecies[rHash] = []
                allSpecies[rHash].append(bimol)
                # trimolecular
                # for spec3 in species:
                #    trimol = [spec1,spec2,spec3]
                #    rHash = self.getExplodedSpecies(trimol)
                #    if(not(rHash in allSpecies)): allSpecies[rHash] = []
                #    allSpecies[rHash].append(trimol)

        # search all branches of a given reactions and store
        allReact = dict()
        for reaction in self.reactions:
            rHash = self.getExplodedSpecies(reaction.reactants)
            if rHash not in allReact: allReact[rHash] = []
            allReact[rHash].append(reaction)

        # loop reactions to find missing and present branches
        for (rHash, reactions) in allReact.items():

            # init data branch structure
            dataBranch = dict()
            dataBranch["reactants"] = reactions[0].reactants
            dataBranch["missingBranches"] = []
            dataBranch["presentBranches"] = []

            # store all reaction branches
            branches = []
            for reaction in reactions:
                branches.append(sorted([x.name for x in reaction.products]))

            # get reactant names (same for all reactions)
            reactNames = sorted([x.name for x in reactions[0].reactants])
            reactNamesNoCR = sorted([x.name for x in reactions[0].reactants if (x.name != "CR")])

            # get reaction hash for reactants only
            rHash = self.getExplodedSpecies(reactions[0].reactants)
            branchFound = []
            # loop on possible species combination for the given reactant hash
            for species in allSpecies[rHash]:
                branch = sorted([x.name for x in species])
                totalCharge = sum([x.charge for x in species])
                maxCharge = max([x.charge for x in species])
                hasElectron = ("E" in [x.name for x in species])
                # skip repulsive products
                if totalCharge == 0 and maxCharge != 0 and not hasElectron: continue
                # skip CR in products
                if "CR" in branch: continue
                # skip X+CR->X
                if branch == reactNamesNoCR: continue
                # skip already-found branches
                if branch in branchFound: continue
                # skip same reactant and products
                if branch == reactNames: continue
                # store branches
                if branch in branches:
                    dataBranch["presentBranches"].append(species)
                else:
                    dataBranch["missingBranches"].append(species)
                branchFound.append(branch)

            # copy data to structure
            missingBranch.append(dataBranch)

        return missingBranch

    # ****************
    def makeHtmlIndex(self):
        fname = "htmls/index.html"

        # open file to write
        fout = open(fname, "w")
        # add header
        fout.write(utils.getFile("header.php"))
        fout.write("<p style=\"font-size:30px\">KROME main menu</p>\n")
        fout.write("<br><br>")
        fout.write("<ul>")
        fout.write("<li><a href=\"indexReactions.html\">Reactions</a></li>")
        fout.write("<li><a href=\"indexSpecies.html\">Species</a></li>")
        fout.write("<li><a href=\"indexGraph.html\">Graphs</a></li>")
        fout.write("<li><a href=\"indexMissingReactions.html\">Missing reactions</a></li>")
        fout.write("<li><a href=\"indexMissingBranches.html\">Missing branches</a></li>")
        fout.write("<li><a href=\"multipleRates.html\">Reactions with multiple rates</a></li>")
        fout.write("</ul>")
        fout.write(utils.getFooter("footer.php"))
        fout.close()

    # ****************
    def makeHtmlGraphIndex(self):

        fname = "htmls/indexGraph.html"

        # open file to write
        fout = open(fname, "w")
        # add header
        fout.write(utils.getFile("header.php"))
        fout.write("<p style=\"font-size:30px\">Graphs</p>\n")
        fout.write("<a href=\"index.html\">back</a><br>\n")
        fout.write("<br><br>\n")

        # loop on all available atoms in the network
        for refAtom in self.getAtoms():

            fnamePNG = "pngs/network_atom_" + str(refAtom) + ".png"
            fnameEPS = "epss/network_atom_" + str(refAtom) + ".eps"
            if not os.path.isfile(fnamePNG): continue

            fout.write("<p style=\"font-size:20px\">Graph for " + refAtom + "</p>\n")
            fout.write("<p>&nbsp;Download <a href=\"../" + fnamePNG + "\" target=\"_blank\">PNG</a>")
            fout.write(" <a href=\"../" + fnameEPS + "\" target=\"_blank\">EPS</a></p><br>\n")
            fout.write("<a target=\"_blank\" href=\"../" + fnamePNG + "\">")
            fout.write("<img src=\"../" + fnamePNG + "\" width=\"700px\" style=\"{max-width:500px;}\"></a>\n")
            fout.write("<br><br>\n")

        # add footer
        fout.write(utils.getFooter("footer.php"))
        fout.close()

    # ****************
    # create reaction list index as html page
    def makeHtmlSpeciesIndex(self):

        fname = "htmls/indexSpecies.html"
        header = "<th><th><th><th><th>"

        # open file to write
        fout = open(fname, "w")
        # add header
        fout.write(utils.getFile("header.php"))
        fout.write("<p style=\"font-size:30px\">Species</p>\n")
        fout.write("<a href=\"index.html\">back</a><br>\n")
        fout.write("<br><br>\n")
        # reaction table
        fout.write("<table width=\"60%\">\n")
        fout.write("<tr>" + header + "\n")
        fout.write("<tr><td>name<td>&Delta;H@0K (kJ/mol)<td>&Delta;H@298.15K (kJ/mol)" \
                   + "<td>&alpha; (&Aring;<sup>3</sup>)<td>xsecs?\n")
        fout.write("<tr>" + header + "\n")
        icount = 0
        # loop on reactions
        for mySpecies in sorted(self.getSpecies(), key=lambda x: x.name):
            bgcolor = ""
            if icount % 2 != 0: bgcolor = utils.getHtmlProperty("tableRowBgcolor")
            enthalpy0 = mySpecies.getEnthalpy(self.thermochemicalData, Tgas=1e-40)
            enthalpy298 = mySpecies.getEnthalpy(self.thermochemicalData)
            polarizability = mySpecies.getPolarizability(self.polarizabilityData)
            if polarizability is not None: polarizability /= 1e-24  # cm3->AA3
            hasPhoto = ""
            if len(mySpecies.phrates) > 0: hasPhoto = "yes"
            fout.write("<tr bgcolor=\"" + bgcolor + "\"><td>&nbsp;" + mySpecies.getHrefName() + "<td>" \
                       + str(enthalpy0) + "<td>" + str(enthalpy298) + "<td>" \
                       + str(polarizability) + "<td>" + hasPhoto + "\n")
            icount += 1
        fout.write("<tr>" + header + "\n")
        fout.write("</table>\n")

        # add footer
        fout.write(utils.getFooter("footer.php"))
        fout.close()

    # ****************
    # create reaction list index as html page
    def makeHtmlReactionIndex(self, myOptions):

        fname = "htmls/indexReactions.html"

        tableHeader = "<tr>" + ("<th>" * 30)

        # open file to write
        fout = open(fname, "w")
        # add header
        fout.write(utils.getFile("header.php"))
        fout.write("<p style=\"font-size:30px\">Reactions</p>\n")
        fout.write("<a href=\"index.html\">back</a><br>\n")

        for variable in myOptions.getRanges().keys():
            fout.write("<a href=\"allRates_" + variable + "_0.html\">All rates with <b>" \
                       + variable + "</b></a><br>\n")
        fout.write("<br><br>\n")
        # reaction table
        fout.write("<table width=\"50%\">\n")
        fout.write(tableHeader + "\n")
        icount = 0
        # reactions sorted by unsorted reaction hash (in the html list)
        reactionsSorted = sorted(self.reactions, key=lambda x: x.getReactionHashUnsorted())
        # loop on reactions
        for myReaction in reactionsSorted:
            bgcolor = ""
            if icount % 2 != 0: bgcolor = utils.getHtmlProperty("tableRowBgcolor")
            fout.write("<tr valign=\"baseline\" bgcolor=\"" + bgcolor + "\">" \
                       + myReaction.getReactionHtmlRow() + "\n")
            icount += 1
        fout.write(tableHeader + "\n")
        fout.write("</table>\n")

        # add footer
        fout.write(utils.getFooter("footer.php"))
        fout.close()

    # ****************
    # create reaction list index as html page
    def makeHtmlMissingBranchesIndex(self, myOptions):

        fname = "htmls/indexMissingBranches.html"

        kJmol2K = 120.274  # kJ/mol->K

        tableHeader = "<tr>" + ("<th>" * 30)

        # open file to write
        fout = open(fname, "w")
        # add header
        fout.write(utils.getFile("header.php"))
        fout.write("<p style=\"font-size:30px\">Missing branches, &Delta;H/K</p>\n")
        fout.write("<a href=\"index.html\">back</a><br><br><br>\n")

        # sort reactions by the name of the first reactant (reactants are also sorted)
        sortedMissingBranch = sorted(self.getMissingBranch(),
                                     key=lambda x: sorted([sp.name for sp in x["reactants"]]))

        fout.write("<table width=\"60%\">\n")
        icount = 0
        for dataBranch in sortedMissingBranch:
            reactants = dataBranch["reactants"]
            reactantsEnthalpy = [x.getEnthalpy(self.thermochemicalData) for x in reactants]
            row = "<td>" + str(icount + 1) + "<td>" \
                  + ("<td>+<td>".join(sorted([x.getHtmlName() for x in reactants])))

            row += ("<td>" * (20 - row.count("<td>")))

            # background color
            bgcolor = utils.getHtmlProperty("tableRowBgcolor")

            # add an arrow at the end of the row
            row += "<td><td>&rarr;"
            # write row to file
            fout.write("<tr valign=\"baseline\" bgcolor=\"" + bgcolor + "\">" + row + "\n")

            # count <td> elements
            ntds = row.count("<td>")

            for listName in ["presentBranches", "missingBranches"]:
                for products in dataBranch[listName]:
                    bgcolor = ""  # this row has no bgcolor
                    productsEnthalpy = [x.getEnthalpy(self.thermochemicalData) for x in products]
                    bgcolor = ""  # this row has no bgcolor
                    # create row from products name
                    rowx = "<td>&rarr;<td>" \
                           + ("<td>+<td>".join(sorted([x.getHtmlName() for x in products])))
                    # add n-1 td as offset
                    row = ("<td>" * (ntds - 1)) + rowx + ("<td>" * (10 - rowx.count("<td>")))
                    status = listName.replace("Branches", "")
                    if status == "missing": status += " &#9888;"

                    tdEnthalpy = ""
                    if None in (productsEnthalpy + reactantsEnthalpy):
                        tdEnthalpy += "<td>missing enthalpy data"
                    else:
                        DeltaH = sum(productsEnthalpy) - sum(reactantsEnthalpy)
                        DeltaH_K = kJmol2K * DeltaH
                        tdEnthalpy += "<td>" + utils.htmlExp(DeltaH_K)
                        if (DeltaH_K < 0e0):
                            tdEnthalpy += "<td>&#10004;&#10004;"
                        elif (DeltaH_K >= 0 and DeltaH_K < 1e4):
                            tdEnthalpy += "<td>&#10004;"
                        else:
                            tdEnthalpy += "<td>&#10006;"

                    fout.write("<tr valign=\"baseline\" bgcolor=\"" + bgcolor + "\">" + row \
                               + "<td style=\"font-size:10px;\">" \
                               + status.upper() + tdEnthalpy + "\n")

            # dataBranch["presentBranches"] = []
            icount += 1
        fout.write("</table>\n")

        # add footer
        fout.write(utils.getFooter("footer.php"))
        fout.close()

    # ****************
    # create reaction list index as html page
    def makeHtmlMissingReactionIndex(self, myOptions):

        fname = "htmls/indexMissingReactions.html"

        kJmol2K = 120.274  # kJ/mol->K

        tableHeader = "<tr>" + ("<th>" * 30)

        # open file to write
        fout = open(fname, "w")
        # add header
        fout.write(utils.getFile("header.php"))
        fout.write("<p style=\"font-size:30px\">Missing reactions*</p>\n")
        fout.write("<p style=\"font-size:10px\">*if reactants are present no check for missing branches!</p>\n")
        fout.write("<a href=\"index.html\">back</a><br>\n")

        rtypes = {1: "unimolecular", 2: "bimolecular", 3: "3-body"}
        for (nreact, rname) in rtypes.items():
            fout.write("<br><br>\n")
            fout.write("<p style=\"font-size:20px\">" + rname.title() + ", &Delta;H/K</p>\n")
            # reaction table
            fout.write("<table width=\"60%\">\n")
            fout.write(tableHeader + "\n")
            icount = 0
            # loop on reactions
            for reaction in self.missingReactions:
                # get reactants
                reactants = reaction["reactants"]
                # consider only number of reactants for the current reaction type
                if len(reactants) != nreact: continue
                # background color
                bgcolor = utils.getHtmlProperty("tableRowBgcolor")
                # create html row
                row = "<td>" + str(icount + 1) + "<td>" \
                      + ("<td>+<td>".join([x.getHtmlName() for x in reactants]))
                # if no branches available (given the network, skip)
                if (len(reaction["products"]) == 0): continue

                # get DH for each reactant
                reactantsEnthalpy = [x.getEnthalpy(self.thermochemicalData) for x in reactants]

                # add an arrow at the end of the row
                row += "<td><td>&rarr;"
                # write row to file
                fout.write("<tr valign=\"baseline\" bgcolor=\"" + bgcolor + "\">" + row + "\n")
                # count <td> elements
                ntds = row.count("<td>")
                # loop on products
                for products in reaction["products"]:
                    productsEnthalpy = [x.getEnthalpy(self.thermochemicalData) for x in products]
                    bgcolor = ""  # this row has no bgcolor
                    # create row from products name
                    rowx = "<td>&rarr;<td>" + ("<td>+<td>".join([x.getHtmlName() for x in products]))
                    # add n-1 td as offset
                    row = ("<td>" * (ntds - 1)) + rowx + ("<td>" * (10 - rowx.count("<td>")))
                    if None in (productsEnthalpy + reactantsEnthalpy):
                        row += "<td>missing enthalpy data"
                    else:
                        DeltaH = sum(productsEnthalpy) - sum(reactantsEnthalpy)
                        DeltaH_K = kJmol2K * DeltaH
                        row += "<td>" + utils.htmlExp(DeltaH_K)
                        if DeltaH_K < 0e0:
                            row += "<td>&#10004;&#10004;"
                        elif DeltaH_K >= 0 and DeltaH_K < 1e4:
                            row += "<td>&#10004;"
                        else:
                            row += "<td>&#10006;"
                    # write html row to file
                    fout.write("<tr valign=\"baseline\" bgcolor=\"" + bgcolor + "\">" + row + "\n")
                icount += 1

            fout.write(tableHeader + "\n")
            fout.write("</table>\n")

        # add footer
        fout.write(utils.getFooter("footer.php"))
        fout.close()

    # *******************
    def makeHtmlAllRates(self, myOptions):

        # loop on variables
        for variable in myOptions.getRanges().keys():

            # default file handler
            fout = None

            # count reactions with variable
            plotsNumber = len([x for x in self.reactions if (x.hasVariable(myOptions, variable))])

            # base of relative path
            fnameBaseRel = "allRates_" + variable + "_"
            # base of "absolute" path
            fnameBase = "htmls/" + fnameBaseRel
            # init counters
            icount = 0
            fileCount = 0
            maxPlotPerPage = 100
            # get pages number
            pagesNumber = int(plotsNumber / maxPlotPerPage) + 1
            # loop on reactions
            for myReaction in self.reactions:
                # skip when the variable is not in the reaction rate
                if not myReaction.hasVariable(myOptions, variable): continue
                if icount % maxPlotPerPage == 0:
                    if fout is not None:
                        fout.write("</table>\n")
                        fout.write(utils.getFooter("footer.php"))
                        fout.close()

                    # prepare a file for each variable
                    fname = fnameBase + str(fileCount) + ".html"

                    # open file to write
                    fout = open(fname, "w")
                    fout.write(utils.getFile("header.php"))
                    fout.write("<p style=\"font-size:30px\">All rate plots with <b>" + variable
                               + "</b> (page " + str(fileCount + 1) + "/" + str(pagesNumber) + ")</p>\n")

                    # multi-page menu back arrow
                    if fileCount - 1 >= 0:
                        pagesMenu = "<a href=\"" + fnameBaseRel + str(fileCount - 1) \
                                    + ".html\">&lt;</a> "
                    else:
                        pagesMenu = "&lt; "
                    # multi-page menu links
                    for ipage in range(pagesNumber):
                        # link only if current page
                        if ipage != fileCount:
                            pagesMenu += "<a href=\"" + fnameBaseRel + str(ipage) + ".html\">" \
                                         + str(ipage + 1) + "</a> "
                        else:
                            pagesMenu += "[" + str(ipage + 1) + "] "
                    # multi-page menu forward arrow
                    if fileCount + 1 < pagesNumber:
                        pagesMenu += "<a href=\"" + fnameBaseRel + str(fileCount + 1) \
                                     + ".html\">&gt;</a>"
                    else:
                        pagesMenu += "&gt;"

                    fout.write(pagesMenu + "<br>")
                    fout.write("<a href=\"indexReactions.html\">back</a><br>\n")
                    fout.write("<br><br>\n")
                    # reaction table
                    fout.write("<table>\n")
                    fileCount += 1

                if icount % 1 == 0: fout.write("<tr><td>")

                fnamePNG = "../pngs/rate_" + str(myReaction.getReactionHash()) + "_" + variable + ".png"
                linkURL = "<a href=\"rate_" + myReaction.getReactionHash() \
                          + ".html\">details</a> for " + myReaction.getVerbatimHtml()
                fout.write("<img src=\"" + fnamePNG + "\" width=\"700px\" alt=\"&#9888; MISSING: " \
                           + myReaction.getVerbatim() + "\"><br>" + linkURL)
                fout.write("<tr height=\"10px\"><td>")
                icount += 1

            if fout is not None:
                fout.write("</table>\n")

                # add footer
                fout.write(utils.getFooter("footer.php"))
                fout.close()

    # *******************
    # prepares HTML documentation for list of rate divided by number of intervals
    def makeHtmlMultipleRates(self, myOptions):

        # divide rates by number of intervals
        multipleRates = dict()
        for reaction in self.reactions:
            # get number of interval
            intervalsNumber = len(reaction.rate)
            # initialize the list for the corresponding dictionary key
            if intervalsNumber not in multipleRates: multipleRates[intervalsNumber] = []
            # add the reaction to dict
            multipleRates[intervalsNumber].append(reaction)

        # prepare a file for each variable
        fname = "htmls/multipleRates.html"

        # standard header for rates
        tableHeader = "<tr>" + ("<th>" * 31)

        # open HTML file to write
        fout = open(fname, "w")
        # add header
        fout.write(utils.getFile("header.php"))
        fout.write("<a href=\"index.html\">back</a><br>\n")
        # loop on multiple rates
        for intervalsNumber in sorted(multipleRates.keys())[::-1]:
            reactionBlock = multipleRates[intervalsNumber]
            fout.write("<br><br>\n")
            # write plural if necessary
            intervalString = "interval" + ("s" if (intervalsNumber > 1) else "")
            # table title
            fout.write("<p style=\"font-size:20px\">Reactions with " + str(intervalsNumber) + " " \
                       + intervalString + "</p>\n")
            # open reaction table
            fout.write("<table width=\"50%\">\n")
            fout.write(tableHeader + "\n")
            # reactions sorted by unsorted reaction hash (in the html list)
            reactionsSorted = sorted(reactionBlock, key=lambda x: x.getReactionHashUnsorted())

            icount = 0
            # loop on reactions
            for myReaction in reactionsSorted:
                bgcolor = ""
                if icount % 2 != 0: bgcolor = utils.getHtmlProperty("tableRowBgcolor")
                fout.write("<tr valign=\"baseline\" bgcolor=\"" + bgcolor + "\">" \
                           + myReaction.getReactionHtmlRow(mode="joints") + "\n")
                icount += 1
            fout.write(tableHeader + "\n")
            fout.write("</table>\n")

        # add footer
        fout.write(utils.getFooter("footer.php"))
        fout.close()

    # ************************
    # create a subnetwork using given options (as subkida.py does)
    def subNetwork(self, myOptions):

        # skip subnetwork if option not present
        if not hasattr(myOptions, 'suboptions'): return

        # expected options divided by type
        listOptions = {"skipAtoms": "", "useAtoms": "", "skipSpecies": "",
                       "skipString": "", "useSpecies": ""}  # comma-separted list
        listPipeOptions = {"skipRateString": ""}  # pipe-separated list
        floatOptions = {"Tmin": -1e99, "Tmax": 1e99, "maxAtoms": 999}
        boolOptions = {"cations": True, "anions": True, "skipTlimitsSingle": False}
        stringOptions = {"outputFile": "subNetwork.dat"}

        # store all keys
        keys = list(listOptions.keys())
        keys += list(listPipeOptions.keys())
        keys += list(floatOptions.keys())
        keys += list(boolOptions.keys())
        keys += list(stringOptions.keys())

        # store options
        fullOptions = dict()
        plainOptions = dict()

        # open options file
        fh = open(myOptions.suboptions)
        for row in fh:
            srow = row.strip()
            if srow == "": continue
            if srow.startswith("#"): continue

            # read options
            (option, value) = [x.strip() for x in srow.split("=")]
            plainOptions[option] = value
            # divide options by type
            if option in listOptions:
                fullOptions[option] = [x.strip() for x in value.split(",") if (x != "")]
            elif option in listPipeOptions:
                fullOptions[option] = [x.strip() for x in value.split("|") if (x != "")]
            elif option in floatOptions:
                fullOptions[option] = float(value)
            elif option in boolOptions:
                sval = value.lower().strip()
                fullOptions[option] = (sval == "t" or sval == "true")
            elif option in stringOptions:
                fullOptions[option] = value
            else:
                print("ERROR: unknown option " + option + " in " + myOptions.suboptions)
                print(" options available are (case sensitive):")
                optionsNames = listOptions.keys() + floatOptions.keys() \
                               + boolOptions.keys() + stringOptions.keys() + listPipeOptions.keys()
                print(", ".join(optionsNames))
                sys.exit()
        fh.close()

        # loop on options keys to fix missing
        for k in keys:
            # replace with empty if option missing
            if k not in fullOptions.keys():
                print("WARNING: " + k + " option is missing in " + myOptions.suboptions + ", replaced with empty")
                fullOptions[k] = ""

        # loop to subselect species
        speciesOK = []
        for species in self.getSpecies():
            # check atoms
            hasSkipAtoms = hasNotUseAtoms = False
            for atom in species.atoms:
                if atom in fullOptions["skipAtoms"]: hasSkipAtoms = True
                if atom not in fullOptions["useAtoms"]: hasNotUseAtoms = True
            if hasNotUseAtoms and fullOptions["useAtoms"] != []: continue
            if hasSkipAtoms: continue

            # check species names
            if species.name in fullOptions["skipSpecies"]: continue

            # include only these species
            if fullOptions["useSpecies"] != []:
                if species.name not in fullOptions["useSpecies"]: continue

            # check special strings
            hasSkipString = False
            for skipString in fullOptions["skipString"]:
                if skipString in species.name: hasSkipString = True
            if hasSkipString: continue

            # check max atoms
            explodedNoSigns = [x for x in species.exploded if (not (x in ["+", "-"]))]
            if len(explodedNoSigns) > fullOptions["maxAtoms"]: continue

            # check ions
            if species.charge > 0 and not (fullOptions["cations"]): continue
            if species.charge < 0 and not (fullOptions["anions"]): continue

            # append OK species
            speciesOK.append(species)

        nameOK = sorted([x.name for x in speciesOK])
        print("species included in the subnetwork (" + str(len(speciesOK)) + "):", (", ".join(nameOK)))

        # print species not found
        if fullOptions["useSpecies"] != []:
            speciesNotFound = [x for x in fullOptions["useSpecies"] if x not in nameOK]
            print("species from " + myOptions.suboptions + " not found:", (", ".join(speciesNotFound)))

        # get all species names
        allNames = [x.name for x in speciesOK]
        reactionsOK = []
        # loop on reactions
        for reaction in self.reactions:
            speciesNames = [x.name for x in (reaction.reactants + reaction.products)]
            # check if all selected species are present
            hasAllSpecies = True
            for species in speciesNames:
                if species not in allNames:
                    hasAllSpecies = False
                    break
            if not hasAllSpecies: continue

            # check special strings in rates
            if fullOptions["skipRateString"] != []:
                # check special strings in rates
                hasSkipRateString = False
                # loop on strings
                for skipString in fullOptions["skipRateString"]:
                    # loop on rates f90 expressions
                    for krate in reaction.rate:
                        # check if string in rate expression
                        if skipString in krate:
                            hasSkipRateString = True
                            break
                # skip reaction if string matches
                if hasSkipRateString: continue

            # create a copy of the reaction to get ranges
            reactionOK = copy.copy(reaction)
            reactionOK.rate = []
            reactionOK.Tmin = []
            reactionOK.Tmax = []

            # replace string in Trange
            TrangeExtras = [x.lower() for x in [".lt.", ".gt.", ".ge.", ".le.", ">", "<", "="]]

            # loop on rate ranges
            for i in range(len(reaction.rate)):
                # check Tmin limit
                if reaction.Tmin[i] is not None:
                    # replace Trange string from Tmin
                    repTmin = reaction.Tmin[i].lower()
                    for TrangeEx in TrangeExtras:
                        repTmin = repTmin.replace(TrangeEx, "")
                    # check Tmin limit
                    if float(repTmin) < fullOptions["Tmin"]: continue
                # check Tmax limit
                if reaction.Tmax[i] is not None:
                    # replace Trange string from Tmin
                    repTmax = reaction.Tmax[i].lower()
                    for TrangeEx in TrangeExtras:
                        repTmax = repTmax.replace(TrangeEx, "")
                    # check Tmax limit
                    if float(repTmax) > fullOptions["Tmax"]: continue
                # append limits
                reactionOK.rate.append(reaction.rate[i])
                reactionOK.Tmin.append(reaction.Tmin[i])
                reactionOK.Tmax.append(reaction.Tmax[i])
            # append reaction to OK list
            reactionsOK.append(reactionOK)

        # write output
        fout = open(fullOptions["outputFile"], "w")
        # header
        fout.write("#########################\n")
        fout.write("#file automatically generated with DOCMAKE on " + utils.getCurrentTime() + "\n")
        fout.write("#changeset: " + utils.getChangeset()[:7] + "\n")
        fout.write("#using the following options\n")
        for (option, value) in plainOptions.items():
            fout.write("# " + option + " = " + value + "\n")
        fout.write("\n\n")

        icount = 0
        # write reactions in KROME format
        for reaction in reactionsOK:
            useTlims = not (fullOptions["skipTlimitsSingle"] and (len(reaction.rate) == 1))
            fout.write(reaction.getKROMEformat(icount, includeTlimits=useTlims))
            icount += len(reaction.rate)

        fout.close()

        # some message
        print("subnetwork saved to " + fullOptions["outputFile"] + " with " + str(icount) + " reactions")
