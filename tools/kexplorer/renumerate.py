import sys
import fileinput

##############################################################################

#renumerate reactions in network file
fname = sys.argv[1].strip()
cnt = 1

for line in fileinput.input(fname,inplace=True):
	#all info for same reaction is separated by empty line
    if "@" in line or "#" in line or line == "\n":
        print line,
    else:
        linesplit = line.strip().split(",")
        linelist = [str(cnt)] + linesplit[1:]
        newline = ",".join(linelist)
        print newline
        cnt = cnt + 1
