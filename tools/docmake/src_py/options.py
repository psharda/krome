class options:

    # *******************
    # constructor
    def __init__(self, fileName):

        optionsData = dict()

        # options that need to be appended (multiple)
        optionAppend = ["range"]

        # initialize list for options to be appended
        for opt in optionAppend:
            optionsData[opt] = []

        # open options file
        fh = open(fileName)
        # loop on file lines
        for row in fh:
            srow = row.strip()
            # skip blank lines and comments
            if srow == "": continue
            if srow.startswith("#"): continue
            # get column-separated options
            (name, value) = [x.strip() for x in srow.split(":") if (x != "")]

            # append options to be appended (optionAppend list)
            if name in optionAppend:
                optionsData[name].append(value)
            else:
                optionsData[name] = value

        # convert option string name into attribute name
        for (k, v) in optionsData.items():
            setattr(self, k, v)

        # read latex options from json file
        if "latexoptions" in optionsData:
            latexoptions.load(optionsData["latexoptions"])

    # *******************
    # get dictionary of ranges
    def getRanges(self):

        # get ranges from option file
        varRanges = dict()
        for rng in self.range:
            # store range name and range limits
            (rangeName, rangeValue) = [x.strip() for x in rng.split("=")]
            varRanges[rangeName] = [float(x) for x in rangeValue.split(",")]

        return varRanges


# *******************
class latexoptions:
    @classmethod
    def load(cls, filename):
        import json

        # read latex options from json file
        with open(filename) as f:
            optionsData = json.load(f)
        for key, value in optionsData.items():
            if type(value) == dict:
                value = {str(k): str(v) for k, v in value.items()}
            elif type(value) == list:
                value = [(str(k), str(v)) for k, v in value]
            setattr(cls, key, value)
