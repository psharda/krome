#random network generator

from random import randint,random
import sys

nrea = 2000 #number of reactions to generate
nspec = 10 #number of random species (names=FK1,FK2,...)
emin = -11. #max rate coefficient, log (kmin=10**emin)
emax = -8. #max rate coefficient, log (kmax=10**emax)
kmax = 1000000 #maximum number of attempts to generate a reacion 
step = 1.5 #decreasing step of species distribution (step=1 equiprobable species)


ints = []
passo = step
if(step<1.): 
	print "ERROR: step must be equal or larger than 1!"
	sys.exit()
mult = 100*passo
for i in range(nspec):
	mult /= passo
	if(mult<1.):
		print "ERROR: mult<1, decrease step"
		sys.exit()
	#print i,int(mult)
	for j in range(int(mult)):
		ints.append(i+1)


fsrea = []
for i in range(nrea):
	for k in range(kmax):
		sp = []
		for j in range(4):
			pick = ints[randint(0,len(ints)-1)]
			sp.append('FK'+str(pick))
			sp = sorted(sp[:2])+sorted(sp[2:])
		fs = "".join(sp)
		if(not(fs in fsrea)):
			break
		if(k>=kmax-1):
			print "max number of attempts reached!",kmax
			sys.exit()
	fsrea.append(fs)
		
	kk = pow(1e1, random()*(emax-emin)+emin)
	print ",".join([str(i+1),sp[0],sp[1],"",sp[2],sp[3],"","","10.","1e4",str(kk)])
