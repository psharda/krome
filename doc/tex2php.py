from latexmath2png import tex2png
import random


#***********************************
class doc():
	sects = []
	#***********************************
	def addsect(self,sect):
		self.sects.append(sect)
#***********************************
class sect():
	subsects = []
	title = "[notitle]"
	number = 0
	text = ""
	eqns = []
	tabs = []
	verbs = []
	#***********************************
	def __init__(self,title,number):
		self.title = title
		self.number = number
		print number,title
	#***********************************
	def addsub(self,sub):
		self.subsects.append(sub)
	#***********************************
	def addeqn(self,eqn):
		self.eqns.append(sub)
	#***********************************
	def addtab(self,tab):
		self.eqns.append(tab)
	#***********************************
	def addverb(self,verb):
		self.verbs.append(verb)
	#***********************************
	def addtext(self,text):
		self.text += text
#*********************************
def makehash():
	return str(int(random.random()*999999)) + str(int(random.random()*999999))

##################################
def clear_dir(folder):
	import os
	for the_file in os.listdir(folder):
	    file_path = os.path.join(folder, the_file)
	    try:
		if os.path.isfile(file_path):
		    os.unlink(file_path)
	    except Exception, e:
		print e

#********************************
def do_eqn(myeqn):
	print "----------------"
	print myeqn
	idx = makehash()
	print idx
	tex2png([myeqn],str(idx),1.4)
	return myeqn.replace(myeqn,"&nbsp;<img src=\"eqn_png/"+idx+"1.png\">&nbsp;")

#***********************************
def make_eqn(txt):
	eqn = False
	myeqn = ""
	for c in txt:
		if(c=="$"):
			if(eqn):
				print "----------------"
				print myeqn
				idx = makehash()
				print idx
				tex2png([myeqn],str(idx),1.4)
				txt = txt.replace(myeqn,"&nbsp;<img src=\"eqn_png/"+idx+"1.png\">&nbsp;")
				eqn = False
			else:
				myeqn = ""
				eqn = True
				continue
		
		if(eqn): myeqn += c
	return txt.replace("$","")
#***********************************
def grep_between(string, cleft, cright):
	store = False
	btw = ""
	for c in string:
		if(cleft==c): 
			store = True
			continue
		if(cright==c): break
		if(store): btw += c
	return btw
#***********************************
#***********************************
#***********************************


fname = "doc.tex"
fh = open(fname, "rb")


clear_dir("eqn_png")

isect = isub = isubsub = 0
position = "pre"
obj = "txt"

php_header = ""
php_footer = ""

toc = php_header + "<ul>"
out = []
for row in fh:
	srow = row.strip()
	slow = srow.lower()
	if("\\end{document}" in slow):
		break
	if("\\begin{document}" in slow):
		position = "doc"
		continue
	if(position=="pre"): continue

	if(len(srow)>0): 
		if(srow[0]=="%"): continue
	else:
		if(isect>0): out[isect-1] += "<br>\n"

	if("\\section" in slow):
		isect += 1
		title = grep_between(slow,"{","}")
		out.append(php_header)
		if(position!="doc"): out[isect-1] += "</p>\n\n"
		fulltitle = str(isect)+" "+title
		#if(isubsub>1): toc += "</ul></li>"
		if(isubsub>1): toc += "</li>\n"
		if(isub>1): toc += "</ul></li>\n\n"
		if(isect>1): toc += "</ul></li>\n\n"
		toc += "<li><a href=\"out"+str(isect)+".html\">"+fulltitle+"</a>\n"
		toc += "<ul>\n"
		out[isect-1] += "<div class=\"section\">"+fulltitle+"</p></div>\n\n"
		out[isect-1] += "<p>"
		isub = 0	
		position = "sect"
		continue
	if(position=="doc"): continue

	if("\\subsection" in slow):
		isub += 1
		title = grep_between(slow,"{","}")
		fulltitle = str(isect)+"."+str(isub)+" "+title
		out[isect-1] += "</p>\n\n"
		out[isect-1] += "<p>"+fulltitle+"</p>\n\n"
		isubsub = 0
		position = "subsect"
		out[isect-1] += "<p>"
		if(isubsub>1): toc += "</li>\n"
		if(isub>1): toc += "</ul></li>\n"
		toc += "<li><a href=\"out"+str(isect)+".html\">"+fulltitle+"</a>\n"
		toc += "<ul>\n"

		continue

	if("\\subsubsection" in slow):
		isubsub += 1
		title = grep_between(slow,"{","}")
		fulltitle = str(isect)+"."+str(isub)+"."+str(isubsub)+" "+title
		out[isect-1] += "</p>\n\n"
		out[isect-1] += "<p>"+fulltitle+"</p>\n\n"
		position = "subsubsect"
		out[isect-1] += "<p>"
		if(isubsub>1): toc += "</li>\n"
		toc += "<li><a href=\"out"+str(isect)+".html\">"+fulltitle+"</a>\n"
		continue

	if("\\subsubsection" in slow):
		isubsub += 1
		title = grep_between(slow,"{","}")
		out[isect-1] += "</p>\n\n"
		out[isect-1] += "<p>"+str(isect)+"."+str(isub)+"."+str(isubsub)+" "+title+"</p>\n\n"
		position = "subsubsect"
		out[isect-1] += "<p>"
		continue

	if("\\begin{itemize}" in slow):
		out[isect-1] += "<ul>\n"
		obj = "list"
		continue

	if("\\end{itemize}" in slow):
		out[isect-1] += "</ul>\n\n"
		obj = "txt"
		continue

	if("\\begin{enumerate}" in slow):
		out[isect-1] += "<ol>\n"
		obj = "list"
		continue

	if("\\end{enumerate}" in slow):
		out[isect-1] += "</ol>\n\n"
		obj = "txt"
		continue

	if("\\item" in slow):
		if(obj=="item"): out[isect-1] += "<br></li>\n"
		out[isect-1] += "<li>\n\n"
		obj = "item"
		srow = srow.replace("\\item","")

	if("\\begin{alltt}" in slow):
		out[isect-1] += "<div style=\"width: 500px; border: 2px solid black\">\n"
		out[isect-1] += "<p style=\"margin: 0; padding: 10px\">\n"
		objold = obj
		obj = "alltt"
		continue

	if("\\end{alltt}" in slow):
		out[isect-1] += "</p></div>\n\n"
		obj = objold
		continue

	if("\\begin{equation}" in slow):
		obj = "eqn"
		out[isect-1] += "<br>"
		continue

	if("\\end{equation}" in slow):
		obj = "txt"
		continue
	
	if(obj=="eqn"): srow = do_eqn(srow)+"<br>"

	srow = make_eqn(srow)
	out[isect-1] += srow
	if(obj=="alltt"): out[isect-1] += "<br>"

i = 0

fout = open("toc.html","w")
fout.write(toc)
fout.close()

for txt in out:
	i += 1
	txt += php_footer
	print "write page ",i
	fout = open("out"+str(i)+".html","w")
	fout.write(txt)
	fout.close()

