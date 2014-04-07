from scipy.optimize import fsolve
from scipy import interpolate
import math,sys

#THIS FILE CONTAINS THE CLASS FOR CREATING COOLING TABLES
#written by Tommaso Grassi and the KROME team (Apr3,2014)

############################################
# PLEASE DO NOT MODIFY UNDER THIS LINE IF
# YOU ARE NOT SURE OF WHAT YOU ARE DOING!
############################################


class krome_tabcool:

	#some constants
	kboltzmann = 1.380658e-16 #erg/K
	hplanck = 6.6260755e-27 #erg*s
	clight = 2.99792458e10 #cm/s

	#define some dictionaries and lists
	energy = dict() #energy in K, key=level number
	g = dict() #multeplicity, key=level number
	A = dict() #Aij coefficient, key=(i,j) tupla
	B = dict() #B coefficient, key=(i,j) tupla
	Cfit = dict() #Collisional coeff, key=collider; this is a dicionary of dictionaries
	colliders = [] #list of colliders

	#########################
	def is_number(self,s):
		try:
			float(s)
			return True
		except ValueError:
			return False

	##############################
	def load_BASECOL(self,fname,collider,verbose=False):
		colliders = self.colliders
		energy = self.energy
		Cfit = self.Cfit
		g = self.g
		A = self.A
		B = self.B
		nu = dict() #frequency in 1/s, key=level number
		
		Cfit[collider] = dict()
		if(not(collider in colliders)): colliders.append(collider)
		if(verbose): print "Reading from "+fname
		#read data from file
		fh = open(fname)
		icount = 0
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			if(not(self.is_number(srow[0]))): continue
			icount += 1
			if(icount==1):
				number_of_temps = int(srow)
				continue

			if(icount==2):
				number_of_trans = int(srow)
				if(verbose): print "Number of transitions:", number_of_trans
				continue
			if(icount==7):
				while("  " in srow):
					srow = srow.replace("  "," ")
				temps = [float(x) for x in srow.split(" ")]
				if(len(temps)!=number_of_temps):
					sys.exit("ERROR: wrong number of temperature values!")
				continue
			if(icount>7):
				while("  " in srow):
					srow = srow.replace("  "," ")
				arow = srow.split(" ")
				level_low = int(arow[0]) - 1
				level_up = int(arow[1]) - 1
				vals = [float(x) for x in arow[4:]]
				if(len(vals)!=number_of_temps):
					sys.exit("ERROR: wrong number of rate values!")
				Cfit[collider][(int(level_up), int(level_low))] = {"T":temps, "val":vals}
				continue

	###########################
	def load_LAMDA(self,fname,verbose=False):
		colliders = self.colliders
		energy = self.energy
		Cfit = self.Cfit
		g = self.g
		A = self.A
		B = self.B
		nu = dict() #frequency in 1/s, key=level number

		if(verbose): print "Reading from "+fname
		#read data from file
		fh = open(fname)
		icount = 0
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			if(srow[0]=="!"): continue
			srow = srow.split("!")[0].strip()
			icount += 1
			#name of the molecule
			if(icount==1): 
				molecule_name = srow
				if(verbose): print "Molecule name:", molecule_name
				continue
			#molecular weight
			if(icount==2):
				molecular_weight = float(srow)
				continue
			#number of levels
			if(icount==3):
				number_of_levels = int(srow)
				if(verbose): print "Number of levels:", number_of_levels
				continue
			#data of the levels
			if(icount>=4 and icount<4+number_of_levels):
				while("  " in srow):
					srow = srow.replace("  "," ")
				arow = srow.split(" ")
				level_number = int(arow[0])-1
				energy[level_number] = float(arow[1]) * 1.42879e0  #cm-1 -> K
				g[level_number] = float(arow[2])
				continue
			#number of radiative transitions
			if(icount==4+number_of_levels):
				number_of_rad_trans = int(srow)
				if(verbose): print "Number of radiative transitions:", number_of_rad_trans
				continue
			#read radiative transition data
			if(icount>4+number_of_levels and icount<=4+number_of_levels+number_of_rad_trans):
				while("  " in srow):
					srow = srow.replace("  "," ")
				(trans, level_up, level_low, level_A, freq, nrg) = srow.split(" ")
				level_up = int(level_up)-1 #upper level
				level_low = int(level_low)-1 #lower level
				level_A = float(level_A)
				freq = float(freq) * 1e9 #Ghz->1/s
				A[(level_up, level_low)] = level_A #store Aij 1/s
				B[(level_up, level_low)] = .5e0 * level_A * self.clight**2 / self.hplanck / (freq)**3
				B[(level_low, level_up)] = B[(level_up, level_low)] * g[level_up] / g[level_low]
				continue
			#number of colliders
			if(icount==4+number_of_levels+number_of_rad_trans+1):
				number_of_colliders = int(srow)
				continue
			#read collider name
			if(icount==4+number_of_levels+number_of_rad_trans+2):
				arow = srow.split(" ")
				#according to LAMDA: collision partner ID and reference.
				# Valid identifications are: 1=H2, 2=para-H2, 3=ortho-H2, 4=electrons, 5=H, 6=He.
				valid_id = ["","H2","H2pa","H2or","e","H","He"]
				partner_id = int(arow[0])
				icount_store = icount
				if(partner_id>len(valid_id)-1 or partner_id==0): 
					sys.exit("ERROR: unknown partner id "+str(partner_id))
				collider = valid_id[partner_id]
				Cfit[collider] = dict()
				if(not(collider in colliders)): colliders.append(collider)
				if(verbose): print "collider fonund:", collider
				continue
			#read number of collisional excitations
			if(icount==4+number_of_levels+number_of_rad_trans+3):
				number_of_coll_trans = int(srow)
				count_colls = 0
	 			continue
			#number of temperatures
			if(icount==4+number_of_levels+number_of_rad_trans+4):
				number_of_temps = int(srow)
	 			continue
			#read temperature values
			if(icount==4+number_of_levels+number_of_rad_trans+5):
				while("  " in srow):
					srow = srow.replace("  "," ")
				temps = [float(x) for x in srow.split(" ")]
				if(len(temps)!=number_of_temps):
					sys.exit("ERROR: wrong number of temperatures "+str(len(temps)))
	 			continue
			#read collisional values
			ntot = 4+number_of_levels+number_of_rad_trans+5
			if(icount>ntot and icount<=ntot+number_of_coll_trans):
				while("  " in srow):
					srow = srow.replace("  "," ")
				arow = srow.split(" ")
				level_up = int(arow[1])-1
				level_low = int(arow[2])-1
				vals = [float(x) for x in arow[3:]]
				if(len(vals)!=number_of_temps):
					sys.exit("ERROR: wrong number of rate values "+str(len(vals)))
				count_colls += 1
				Cfit[collider][(int(level_up), int(level_low))] = {"T":temps, "val":vals}
				#rewind the number of lines to read all colliders
				if(count_colls==number_of_coll_trans):
					if(number_of_colliders>1):
						number_of_colliders -= 1
						icount = icount_store-1
				continue




	###########################
	#load data from krome file with flin
	def load_krome(self,fname,verbose=False):
		colliders = self.colliders
		energy = self.energy
		Cfit = self.Cfit
		g = self.g
		A = self.A
		B = self.B
		nu = dict() #frequency in 1/s, key=level number

		if(verbose): print "Reading from "+fname
		#read data from file
		fh = open(fname)

		#loop on file
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			if(srow[0]=="#"): continue
			srow = srow.split("#")[0]
			srow = srow.split("//")[0]

			#read and store level data: number, energy, multeplicity
			#compute also frequency
			if("level" in srow): 
				arow = srow.split(":")
				level_number = int(arow[0].replace("level","").strip()) #level number (zero-based)
				(level_energy, level_g) = [x.strip() for x in arow[1].split(",")]
				energy[level_number] = float(level_energy) #energy in K
				nu[level_number] = float(level_energy) * self.kboltzmann / self.hplanck #frequency in 1/s
				g[level_number] = float(level_g) #level multeplicity

			#read and store spontaneus transitions: upper level, lower level, Aij coefficient
			if(len(srow.split(","))==3):
				(level_up, level_low, level_A) = [x.strip() for x in srow.split(",")]
				level_up = int(level_up) #upper level
				level_low = int(level_low) #lower level
				level_A = float(level_A) #Aij coefficient 1/s
				A[(level_up, level_low)] = level_A #store Aij
				#compute Bji and Bij from Aij
				B[(level_up, level_low)] = .5e0 * level_A * self.clight**2 / self.hplanck / (nu[level_up] - nu[level_low])**3
				B[(level_low, level_up)] = B[(level_up, level_low)] * g[level_up] / g[level_low]

			#retrive collisional data from flin F90-style function
			# data are: collider name, upper level, lower level, rate in F90-style
			if("flin" in srow):
				(collider, collider_up, collider_down, rate) = srow.split(",",3)[:4]
				#explode flin to obtain Trange and corresponding function values
				arate = rate.replace("d","e").replace("(","").replace(")","").replace(", Tgas","").replace("flin","").split("/, /")
				arate = [x.replace("/","").strip() for x in arate]
				Trange = [float(x) for x in arate[0].split(",")] #store Trange values in a list
				vrate = [float(x) for x in arate[1].split(",")] #store coefficient values in a list
				if(len(Trange)!=len(vrate)): sys.exit("ERROR: Trange must be equal to vrate!")
				#store collider names and define the dictionary of Cij for the given collider
				if(not(collider in colliders)):
					colliders.append(collider)
					Cfit[collider] = dict()
				#store Trange and vrate in the dictionary
				Cfit[collider][(int(collider_up), int(collider_down))] = {"T":Trange, "val":vrate}

		if(verbose):
			print "lines found:",len(energy)
			print "colliders and collisional transitions found:"
			for coll in colliders:
				print " ",coll,len(Cfit[coll])

			print "spontaneus transitions found:",len(A)
			print "Initialization: OK!"

		self.A = A
		self.B = B
		self.Cfit = Cfit
		self.colliders = colliders
		self.energy = energy
		self.g = g
		

	#######################
	#function to compute tau according to NK93
	def ftau(self,lBji,lBij,dv,xj,xi):
		return self.hplanck*self.clight*0.25/math.pi/dv*(xj*lBji-xi*lBij)

	#######################
	#function to beta according to NK93
	def fbeta(self,lBji,lBij,dv,xj,xi):
		mytau = self.ftau(lBji,lBij,dv,xj,xi)
		return 1e0/(1e0+3e0*mytau)

	######################
	#define system of equations as in NK93
	def eqs(self,x,args):
		#read data from additional arguments
		dvdz = args["dvdz"]
		lA = args["A"]
		lB = args["B"]
		lC = args["C"]
		lcolliders = args["colliders"]
		lxcoll = args["xcoll"]
		nlev = args["nlev"]
		Tbb = args["Tbb"]
		energy = args["energy"]
		eq = [0e0 for i in range(nlev)] #initialize equations
		#loop on levels to prepare equations
		for i in range(nlev):
			#prepare first term (j->i)
			p1 = 0e0
			#loop on levels
			for j in range(nlev):
				if(i==j): continue #i=j is not a transition!

				AA1 = 0e0
				#use only available transitions
				if((j,i) in lA):
					AA1 = lA[(j,i)]
					AA1 *= self.fbeta(lB[(i,j)],lB[(j,i)],dvdz,x[i],x[j]) #add shielding

				CC = 0e0
				#sum up (de)excitations from colliders
				for coll in lcolliders:
					if((j,i) in lC[coll]): CC += lC[coll][(j,i)] * lxcoll[coll]

				BI = 0e0
				#stimulated emmission/absorption for when blackbody radiation with Tbb
				if(((j,i) in lA) and (Tbb>0e0)):
					dE = abs(energy[i]-energy[j]) #K
					BI = AA1/(math.exp(dE/Tbb)-1e0)

				#put everything together
				p1 += x[j] * (AA1 + BI + CC)

			#prepare second term (i->j), see comments above
			p2 = 0e0
			for j in range(nlev):
				if(i==j): continue #i=j is not a transition!
				AA2 = 0e0
				if((i,j) in lA):
					AA2 = lA[(i,j)]
					AA2 *= self.fbeta(lB[(j,i)],lB[(i,j)],dvdz,x[j],x[i])

				CC = 0e0
				for coll in lcolliders:
					if((i,j) in lC[coll]): CC += lC[coll][(i,j)] * lxcoll[coll]

				BI = 0e0
				if(((i,j) in lA) and Tbb>0e0):
					dE = abs(energy[i]-energy[j]) #K
					BI = AA2/(math.exp(dE/Tbb)-1e0)

				p2 += (AA2 + BI + CC)

			#build equation for the ith level
			eq[i] = p1 - x[i] * p2
	
		#first equation is replaced by continuity
		eq[0] = sum(x) - 1e0
		return eq

	########################
	#recap the number of lines and transitions
	def recap(self):
		colliders = self.colliders
		Cfit = self.Cfit
		Aij = self.A
		print "*********************"
		print "      DATA RECAP"
		for collider in colliders:
			mytrans = []
			maxlevC = 0
			minlevC = 9999
			for k,v in Cfit[collider].iteritems():
				if(not(k in mytrans)): 
					mytrans.append(k)
					maxlevC = max(list(k)+[maxlevC])
					minlevC = min(list(k)+[minlevC])
			print "*********************"
			print "collider:",collider
			nlev = maxlevC-minlevC+1
			print "collisional transitions (Cij):",len(mytrans)," of ",nlev*(nlev-1)/2,"allowed"
			print "upper level:",maxlevC
			print "lower level:",minlevC
			print
		maxlevA = 0
		minlevA = 9999
		for k,v in Aij.iteritems():
			maxlevA = max(list(k)+[maxlevA])
			minlevA = min(list(k)+[minlevA])
		print "*********************"
		print "spontaneus transitions (Aij):",len(Aij)
		print "upper level:",maxlevA
		print "lower level:",minlevA

		print


			
	#######################
	#compute cooling
	def get_cool(self,Tgas,xcoll,absdvdz=1e99,Tbb=0e0):
		colliders = self.colliders
		energy = self.energy
		Cfit = self.Cfit
		g = self.g
		A = self.A
		B = self.B
		C = dict() #init collision dictionary
		#loop on colliders
		for coll in colliders:
			C[coll] = dict()
			#loop on transitions for the given collider
			for k,v in Cfit[coll].iteritems():
				#interpolate Cij for Tgas
				f = interpolate.interp1d(Cfit[coll][k]["T"], Cfit[coll][k]["val"])
				Tgasi = max(Tgas,Cfit[coll][k]["T"][0])
				Tgasi = min(Tgasi,Cfit[coll][k]["T"][len(Cfit[coll][k]["T"])-1])
				C[coll][k] = f(Tgasi) #store Cij
				deltaE = energy[k[0]]-energy[k[1]] #already K
				#compute and store Cji from Cij
				C[coll][(k[1],k[0])] = C[coll][k] * g[k[0]]/g[k[1]] * math.exp(-deltaE/Tgas)


		#prepare additional arguments
		myargs = {"A":A, "B":B, "C":C,"nlev":len(energy),"dvdz":absdvdz,"colliders":colliders,"xcoll":xcoll,"Tbb":Tbb,"energy":energy}

		#initial guess
		x0 = [1e0/len(energy) for i in range(len(energy))]
	
		#solve the (non)linear system
		res = fsolve(self.eqs, x0,args=myargs,full_output=True)
		if(res[2]!=1): sys.exit("ERROR: "+res[3]) #check for error
		xx = res[0] #copy solution
		#compute cooling as x(j)*Aji*betaji (erg/s/cm3)
		cool = 0e0
		for k,v in A.iteritems():
			kup = k[0]
			klow = k[1]
			beta = self.fbeta(B[(klow,kup)],B[(kup,klow)],myargs["dvdz"],xx[klow],xx[kup])
			cool += A[k] * xx[kup] * (energy[kup] - energy[klow]) * self.kboltzmann
		return cool





