from random import randint,random
import sys

nrea = 2000 #numero di reazioni da generare
nspec = 10 #numero di specie FK1,FK2,...
emin = -11. #esponenziale minimo rate coefficient (10**emin)
emax = -8. #esponenziale massimo rate coefficient (10**emax)
kmax = 1000000 #numero massimo di tentativi per generare una reazione 


ints = []
passo = 1.5
mult = 100*passo
for i in range(nspec):
	mult /= passo
	if(mult<1.):
		print "ERRORE: mult<1, dimunusici passo"
		sys.exit()
	print i,int(mult)
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
			print "numero di tentativi superato!",kmax
			sys.exit()
	fsrea.append(fs)
		
	kk = pow(1e1, random()*(emax-emin)+emin)
	print ",".join([str(i+1),sp[0],sp[1],"",sp[2],sp[3],"","","10.","1e4",str(kk)])
