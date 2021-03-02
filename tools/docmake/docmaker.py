#!/usr/bin/env python
import sys

# move to the python source folder
sys.path.insert(0, "./src_py/")

import network
import options

# check argv
if len(sys.argv) != 2:
    print("USAGE is: python docmaker.py OPTION_FILE")
    sys.exit()

# get filename
fname = sys.argv[1]

# load options
myOptions = options.options(fname)

# load network
myNetwork = network.network(myOptions)

# prepare full documentation
myNetwork.makeAllDoc()
