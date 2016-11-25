#this class compute optically thin photorates for atoms given a radiation flux
# see usage below
from math import sqrt,pi,exp

class photorates:

	#xsecs data
	data = dict()
	#flux data
	dataFlux = []

	#constants
	hplanck = 4.135667662e-15 #eV*s
	kboltzmann = 8.6173303e-5 #eV/K
	clight = 2.9979245800e10 #cm/s
	stefbol = 3.5391853e7 #eV/cm2/s/K4

	#default Draine flux limits
	eminDraine = 6e0 #eV
	emaxDraine = 13.6e0 #eV

	#***************************
	#constructor and load data from Verner
	def __init__(self,vernerFile="photo.dat"):

		#atoms list (atomic number)
		atoms = {1:"H", 2:"He", 6:"C", 7:"N", 8:"O", 14:"Si"}
		#fit parameters
		ckey = ["E_th","E_max","E_0","sigma_0","y_a","P","y_w","y_0","y_1"]
		#load from file
		fh = open(vernerFile,"rb")
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			if(srow.startswith("#")): continue
			#read line
			arow = [x for x in row.split(" ") if(x!="")]
			#dump coefficient into dictionary
			coeff = {ckey[i]:float(arow[2+i]) for i in range(len(ckey))}
			#get atomic number, and electrons number
			(Z,N) = [int(x) for x in arow[:2]]
			#skip missing atoms
			if(not(Z in atoms)): continue
			#create atom name, e.g. C++
			atomName = atoms[Z]+("+"*(Z-N))
			#save fit coefficients into dictionary
			self.data[atomName] = coeff
		fh.close()

	#***************************
	#get xsec for energy/eV, in cm2
	def getXsec(self,atomName,energy):
		data = self.data[atomName]
		#convert coefficient names into variables
		for (k,v) in data.iteritems():
			exec(k + " = "+str(v)+"")
		#check energy threshold
		if(energy<=E_th): return 0e0
		#see Verner+1996, ApJ, 465, 487
		xx = energy/E_0 - y_0
		yy = sqrt(xx**2+y_1**2)
		Fy = ((xx-1.)**2+y_w**2)*yy**(0.5*P-5.5)*(1.+sqrt(yy/y_a))**(-P)

		#Mb->cm2
		return 1e-18*sigma_0*Fy

	#***************************
	#load flux for energy/eV from file, in cm2
	def loadFlux(self,fileName,Emin=0e0,Emax=1e99):
		self.dataFlux = []
		fh = open(fileName,"rb")
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			if(srow.startswith("#")): continue
			(energy,flux) = [float(x) for x in row.split(" ") if(x!="")]
			if(energy<Emin or energy>Emax): flux = 0e0
			self.dataFlux.append([energy,flux])

		fh.close()

	#*****************
	#set Draine flux
	def setDraine(self,nbins=100,emin=None,emax=None):

		#when MIN limit is not set uses default, otherwise overrides
		if(emin==None):
			emin = self.eminDraine
		else:
			self.eminDraine = emin

		#when MAX limit is not set uses default, otherwise overrides
		if(emax==None):
			emax = self.emaxDraine
		else:
			self.emaxDraine = emax

		self.dataFlux = []
		#loop to store flux
		for i in range(nbins):
			xe = i*(emax-emin)/(nbins-1)+emin
			flux = (1.658e6*xe - 2.152e5*xe**2 + 6.919e3*xe**3) * xe * self.hplanck
			self.dataFlux.append([xe,flux])

	#***************
	#integral function of Draine flux
	def integralDraine(self,xe):
		return (1.658e6*xe**3/3. - 2.152e5*xe**4/4. + 6.919e3*xe**5/5.) * self.hplanck

	#***************
	#compute Draine flux
	def getSumDraine(self):
		emin = self.eminDraine
		emax = self.emaxDraine
		return self.integralDraine(emax) - self.integralDraine(emin)

	#***************
	#black body @ Tbb/K, energy/eV, returns eV/cm2/sr
	def fBB(self,Tbb,energy):
			xe = energy
			xexp = xe/self.kboltzmann/Tbb
			#avoid crazy exponential
			if(xexp<3e2 and xe>1e-10):
				return 2e0*xe**3 / self.hplanck**2 / self.clight**2 / (exp(xexp)-1e0)
			return 0e0

	#***************
	#set black body radiation @ Tbb, if no limits uses Draine limits
	def setBB(self,Tbb,Emin=None,Emax=None,nbins=100,verbose=True,normalizeToDraine=True):

		#uses Draine limits if not defined
		if(Emin==None): Emin = self.eminDraine
		if(Emax==None): Emax = self.emaxDraine

		sumBB = 0e0
		#energy bin size (assumes equally-spaced)
		de = (Emax-Emin)/nbins
		self.dataFlux = []
		#loop on bins
		for i in range(nbins):
			#energy
			xe = i*(Emax-Emin)/(nbins-1)+Emin
			#compute BB flux
			flux = self.fBB(Tbb,xe)
			#approximated integral
			sumBB += flux*de
			#store data
			self.dataFlux.append([xe,flux])

		#compute theoretical integral
		sumBB_th = self.stefbol*Tbb**4*self.hplanck
		#include geometry factor
		sumBB *= pi
		fscale = self.getSumDraine()/sumBB*4.
		#print a report when verbose mode on (default)
		if(verbose):
			print "************************"
			print "Tbb (K)",Tbb
			print "BB limits (eV)",Emin,Emax
			print str(sumBB/sumBB_th*1e2)+"% of BB integral"
			print "BB normalization to Draine flux " + str(fscale)
			print "************************"

		#normalize to Draine flux (same integral)
		if(normalizeToDraine):
			self.scaleFlux(fscale)


	#****************
	#scale flux
	def scaleFlux(self,factor):
		self.dataFlux = [[x[0],x[1]*factor] for x in self.dataFlux]

	#*****************
	#compute rate as integral 4pi*J(E)*sigma(E)/E/h dE, 1/s
	def calcRate(self,atomName):

		krate = 0e0
		for i in range(len(self.dataFlux)-1):
			(energyL,fluxL) = self.dataFlux[i]
			(energyR,fluxR) = self.dataFlux[i+1]
			fL = self.getXsec(atomName,energyL)*fluxL/energyL
			fR = self.getXsec(atomName,energyR)*fluxR/energyR
			dE = energyR-energyL
			krate += dE*(fR+fL)/2e0
		return 4.*pi*krate/self.hplanck



