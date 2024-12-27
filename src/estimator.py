#!/usr/bin/env python3
import models
import phymcmc
import numpy
import scipy
import sys


# Base parameter values
#	for fixed and estimated params
base_params_dict = dict(
	# Modelling choices
	# MM = SM or MFM, Likelihood = SSR or ED, [V] = SIN or IV
	SM = False, LikeED = True, IV = True,
	############# Parameters to be estimated (initial guesses)
    # Production rate of V (IV or SIN)
    rho = 10**1.1, # IV/h or SIN/h
	# Production rate of Vrna
	rho_rna = 10**3.1, # vRNA/h
	# Duration of Eclipse, Infectious phases
	tauE = 7., tauI = 27., # h
	################## Parameters held fixed ##################################
	# ---- Infection parameters
	# Number of Eclipse, Infectious compartments (Erlang shape parameter)
	nE = 60, nI = 60,
	# Decay rate of V, Vrna
	c = 0.0573, crna = 0.001, # per h
	# Number of cells in SC, MC infection assays
	Ncells = 1.9e6, # cells
	# Volume of supernatant in SC, MC infection assays
	S = 10., # ml
	# MOI(-1h) in SC to get V(-1), V(0) in SC post-rinse, V(0) in MC
	MOIsc = 3., V0pr= 10.**7.01, V0mc = 10.**2.80, # SIN
	# Vrna(-1) in SC, Vrna(0) in SC post-rinse, Vrna(0) in MC
	Vm1scrna = 10.**11.0, V0prrna = 10.**7.77, V0mcrna = 10.**5.04, # vRNA
	# At each sample time, 0.5 ml or 10 ml is replaced w fresh medium
	frinse = 1.-0.5/10.,
	# --- Std dev of log10(Vsin) across replicates
	# Computed as stdev of list [ dat(t)-mean(dat(t)) ] over all t
	sigma_SC = 0.250, sigma_MC = 0.244, # std dev of log10(Vsin)
	sigma_SC_RNA = 0.229, sigma_MC_RNA = 0.258, # std dev of log10(Vrna)
	# Lower limit of detection of the ED assay in SIN/ml given Vinoc=0.05 ml
	# 	4 replicate/dilution, and 8 serial 10-fold dilutions of 10^[-1,-8]
	# see Quirouette 2025 Methods Eqn (20) FIXME
	LLoD = numpy.log(2.) / (0.05*4*numpy.logspace(-1,-8,8).sum()), # =10^1.494 SIN/ml
	# ---- ED assay parameters
	# Number of cells in each well of the ED assay
	NcellsED = 1.e5, # cells
	# Total volume (virus sample+dilution) in each ED assay well
	Vinoc = 0.05, # ml
	# Random number seed (req'd if SM=True)
	RNseed = 1024, # integer [0,2^32]
)


def get_V0_from_MOI(MOI,tincub,pdic):
	""" Determine V(-1h) in SC based on params so that MOI=3. """
	def diff_actual_expected(V0):
		def ode(t,y):
			T,V = y
			bTV = pdic['beta']/pdic['S']*T*V
			dT = -pdic['gamma']*bTV
			dV = -pdic['c']*V-bTV
			return numpy.hstack((dT,dV))
		# Solve ODE model for tincub (duration of incubation time in h)
		sol = scipy.integrate.solve_ivp(ode,[0.,tincub],[pdic['Ncells'],V0],method='BDF')
		# Diff btwn frac cell uninfected, T(0h), minus expected based on MOI
		return sol.y[0][-1]/pdic['Ncells'] - numpy.exp(-MOI)
	# Guess V0 = [c+b/S*Ncells] * MOI /  [ 1 - exp{-(c+b/S*Ncells)*1h} ] / [g*b/S]  
	V0guess = pdic['c'] + pdic['beta']/pdic['S']*pdic['Ncells']
	V0guess *= MOI / -numpy.expm1(-V0guess*tincub)
	V0guess /= pdic['gamma']*pdic['beta']/pdic['S']
	res = scipy.optimize.root_scalar(diff_actual_expected,method='secant',x0=V0guess,rtol=1e-6)
	if not res.converged:
		diff = abs(diff_actual_expected(res.root))
		print('** WARNING! T/Ncells-exp(-MOI) = %g for:\npvec = %s'%(diff,repr(pdic)),file=sys.stderr)
	return res.root


def get_PVestED(pdic):
	""" Compute P_V->Establishment in ED assay based on parameters in pdic. """
	# If [V]=SIN then P_V->Establishment = 1 (establishment is certain by definition)
	if not pdic['IV']:
		return 1.
	# Equation (5) in Quirouette 2023 (with params for the ED assay)
	# P_V->I = gamma / [ 1 + ( c / {beta*NcellsED/Vinoc} ) ]
	P_V2I = pdic['gamma'] / ( 1. + pdic['c']/(pdic['beta']*pdic['NcellsED']/pdic['Vinoc']) )
	# Equation (30) in Quirouette 2023
	# 0 = 1 - PVest/P_V->I - { 1 / [p*tauI/nI PVest + 1] }^{nI}
	f = lambda pvest: 1. - pvest/P_V2I - (pdic['rho']*pdic['tauI']/pdic['nI']*pvest+1.)**-pdic['nI']
	# Initial guess is that P_V->Establishment = P_V->I
	#	most accurate when burst size (p*tauI) is large
	res = scipy.optimize.root_scalar(f, method='secant', x0=P_V2I, rtol=1e-6)
	if not res.converged:
		diff = abs(f(res.root))
		print('** WARNING! f(PVest) = %g for:\npvec = %s'%(diff,repr(pdic)),file=sys.stderr)
	return res.root


class params(phymcmc.ParamStruct):

	def validate(self):
		""" Assess validity of model parameters. """
		# No estimated parameter is allowed to be negative
		if min(self.vector) < 0:
			raise ValueError('Negative params',repr(self.pardict))

		# Before testing anything, pre-compute quantites we need
		if self.pardict['IV']: # if [V]=IV
			self.pardict['gamma'] = numpy.sqrt(self.pardict['aiv']*self.pardict['biv'])
			self.pardict['rho'] = numpy.sqrt(self.pardict['aiv']/self.pardict['biv'])
		else: # if [V]=SIN
			self.pardict['gamma'] = numpy.sqrt(self.pardict['asin']*self.pardict['bsin'])
			self.pardict['beta'] = numpy.sqrt(self.pardict['asin']/self.pardict['bsin'])
		self.spares = dict(
			kE = 1.*self.pardict['nE']/self.pardict['tauE'],
			dI = 1.*self.pardict['nI']/self.pardict['tauI'],
			# PVestablishment in ED assay, NOT in SC/MC infections
			PVestED = get_PVestED(self.pardict),
			# 1.0 = SC incubation time is 1h (-1h to 0h)
			Vm1sc = get_V0_from_MOI(self.pardict['MOIsc'],1.0,self.pardict),
		)

		# If MM is stochastic model
		if self.pardict['SM'] and self.pardict['RNseed'] >= 2**32:
			raise ValueError('RNseed',repr(self.pardict))
		else:
			# Seed the random number generator with integer
			self.rng = numpy.random.default_rng(round(self.pardict['RNseed']))

		# If V in units of IV
		if self.pardict['IV']:
			if self.pardict['gamma'] > 1.0:
				raise ValueError('gamma>1',repr(self.pardict))
			if self.pardict['rho'] > self.pardict['rho_rna']:
				raise ValueError('rho>rho_rna',repr(self.pardict))
			if self.spares['PVestED'] < (self.pardict['V0pr']/self.pardict['V0prrna']):
				raise ValueError('PVEstED=%g such that V>Vrna in SCpr'%self.spares['PVestED'],repr(self.pardict))
			if self.spares['Vm1sc'] > self.pardict['Vm1scrna']:
				raise ValueError(r'Vsin(-1)>Vrna(-1) in SC', repr(self.pardict))
		# If V in units of SIN
		else:
			if self.pardict['beta'] < 1.e-10:
				raise ValueError(repr(self.pardict))



class estimator(phymcmc.base_model):
	def __init__(self, dat, params):
		self.dat = dat
		self.params = params

	def get_solution(self):
		# Get solution for SC infection
		resSC = models.infection_simulator(self.params,self.dat['tSC'],'SC')
		if resSC is None: # Viv>Vrna
			return None
		# Get solution for MC infection
		resMC = models.infection_simulator(self.params,self.dat['tMC'],'MC')
		if resMC is None: # Viv>Vrna
			return None
		return resSC[1], resSC[2], resMC[1], resMC[2]

	def ln_prior(self,pdic):
		# ln(prior) = ln[ 1 / (beta * gamma * p * prna) ] = -ln[ prod pars ]
		# If we want prior(rho,gamma) = 1/[rho x gamma], we need joint
		#	prior for a,b of prior(a,b) = 1/[a x b]
		if pdic['IV']:
			return -numpy.log(pdic['aiv']*pdic['biv']*pdic['beta']*pdic['rho_rna'])
		return -numpy.log(pdic['asin']*pdic['bsin']*pdic['rho']*pdic['rho_rna'])

	def ln_probability(self,pvec):
		# Set parameter dictionary to new pvec guess
		try:
			self.params.vector = pvec
		except ValueError:
			return float('-inf')

		res = self.get_solution()
		# This is to handle case where Viv >= Vrna
		#	which causes the solver to return None
		if res is None:
			return float('-inf')
		# Otherwise unpack the solution
		scV, scVrna, mcV, mcVrna = res

		# Short-hand for param dictionary
		pdic = self.params.pardict

		# lnprob for Vsin or Viv based on ED or SSR likelihood
		if pdic['LikeED']:
			lnprobV = self.ln_LikeED(scV,mcV)
		else:
			lnprobV = self.ln_LikeSSR(scV,mcV)

		# Short-hand for data
		dat = self.dat

		# lnprob for Vrna in SC (see data_packer)
		#	lnprob = [ log(Vnra_model)[SC_RNA_i] - SC_vRNA ]**2.0 * mul
		#   SC_vRNA = log(vRNA_expt * S) # vRNA
		#   SC_vRNA_norm = -1 / [ 2 * ln(10) * sigma^2 ]
		lnprobVrna = numpy.sum( (numpy.log(scVrna)[dat['SC_RNA_i']]-dat['SC_lnvRNA'])**2. )*dat['SC_lnvRNA_norm']
		# lnprob for Vrna in MC
		lnprobVrna += numpy.sum( (numpy.log(mcVrna)[dat['MC_RNA_i']]-dat['MC_lnvRNA'])**2. )*dat['MC_lnvRNA_norm']

		# lnPrior
		lnprior = self.ln_prior(pdic)

		return lnprobV + lnprobVrna + lnprior


	def ln_LikeSSR(self,scV,mcV):
		pdic = self.params.pardict
		dat = self.dat
		# SC Vsin
		# Convert Cvir from IV->SIN or SIN->SIN (PVext=0 when [V]=SIN)
		lnV = numpy.log( scV * self.params.spares['PVestED'] )
		lnprobV = numpy.sum( (lnV-dat['SC_lnVsin'])**2. )*dat['SC_lnVsin_norm']
		# MC Vsin
		# Convert V from IV->SIN or SIN->SIN (PVext=0 when [V]=SIN)
		lnV = numpy.log( mcV * self.params.spares['PVestED'] )
		# Array MC_lnVsin excludes values below LLoD
		lnprobV += numpy.sum( (lnV-dat['MC_lnVsin'])**2. )*dat['MC_lnVsin_norm']
		# Array MC_lnVsin_nan includes only values below LLoD
		#	only count positive diffs where Vmodel > LLoD
		lnprobV += numpy.sum( numpy.maximum(lnV-dat['MC_lnVsin_nan'],0.)**2. )*dat['MC_lnVsin_norm']
		return lnprobV


	def ln_LikeED(self,scV,mcV):
		dat = self.dat
		lnprobV = 0.
		for expt,V in zip(('SC','MC'),(scV,mcV)):
			# Repeat each Vsin x 3 samples
			Vsin = numpy.repeat(V * self.params.spares['PVestED'], 3)
			# {expt}_sumDnmk = sum_{t,rep} -Vinoc/Vsin * sum_{col} [ Dcol * (ncol-kcol) ]
			lnprobV += numpy.sum(Vsin*dat[f'{expt}_sumVDnmk'])
			# Repeat 3xVsin x 8 columns and reshape
			Vsin = numpy.reshape(Vsin.repeat(8),(-1,8))
			lnprobV += numpy.sum(dat[f'{expt}_k'] * numpy.log1p(-numpy.exp(Vsin*dat[f'{expt}_VD'])))
		return lnprobV



def adjust_pdic_pfit(pdic,pfit):
	# Replace correlated parameters, as needed in pdic
	beta, gamma, rho = pdic.pop('beta'), pdic.pop('gamma'), pdic.pop('rho')
	if pdic['IV']:
		pdic['aiv'] = gamma*rho
		pdic['biv'] = gamma/rho
		pdic['beta'] = beta
		pfit += ['aiv','biv','beta']
	else:
		pdic['asin'] = gamma*beta
		pdic['bsin'] = gamma/beta
		pdic['rho'] = rho
		pfit += ['asin','bsin','rho']
	return pdic, pfit



def data_packer(pdic,data_dir):
	""" Organizes the data so Likelihood computation is fast/efficient. """
	dat = dict()

	# Vrna in SC
	# Discard data at t=-1h,0h because used to set Vrna(-1) & postrinse
	# 	so start from t=1h (skip first 4 rows)
	tmp = numpy.loadtxt(data_dir+'/sH1N1_SC_RNA.dat')[4:,:].T
	# Set sampling intervals for SC
	dat['tSC'] = numpy.unique(numpy.hstack([0.,tmp[0]]))
	# Find indices in tSC when Vrna samples were taken (t=0 is omitted in Vsin)
	dat['SC_RNA_i'] = numpy.hstack([numpy.flatnonzero(dat['tSC']==ts) for ts in tmp[0]])-1
	# Take ln of Vrna*Vol to get ln(vRNA) count rather than concentration
	dat['SC_lnvRNA'] = numpy.log( tmp[1]*pdic['S'] )
	# Denominator of SSR (-2 sigma^2)
	#	since we'll take ln(data) rather than log10(data) for computational
	#	efficiency, sigma_SC_RNA quoted in paper is multiplied by ln(10)
	dat['SC_lnvRNA_norm'] = -0.5 / (numpy.log(10.) * pdic['sigma_SC_RNA'])**2.
	# We will compute lnprob for V_rna in SC as
	#	lnprob = [ log(Vnra_model)[SC_RNA_i] - SC_lnvRNA ]**2.0 * SC_lvRNA_norm

	# Vrna in MC (see SC for comments)
	tmp = numpy.loadtxt(data_dir+'/sH1N1_MC_RNA.dat').T
	dat['tMC'] = numpy.unique(tmp[0])
	dat['MC_RNA_i'] = numpy.hstack([numpy.flatnonzero(dat['tMC']==ts) for ts in tmp[0]])
	dat['MC_lnvRNA'] = numpy.log( tmp[1]*pdic['S'] )
	dat['MC_lnvRNA_norm'] = -0.5 / (numpy.log(10.) * pdic['sigma_MC_RNA'])**2.

	if pdic['LikeED']:
		# SC objects: dilution (D), positive well count (k), total well count (n)
		# Starts at time t=1h (remove t=0 time points)
		SCtmp = numpy.loadtxt(data_dir+'/parsed-SC-data.dat')[3:]
		MCtmp = numpy.loadtxt(data_dir+'/parsed-MC-data.dat')
		# columns are ED assay columns, rows are 3 replicates @ all time point
		for expt,tmp in zip(('SC','MC'),(SCtmp,MCtmp)):
			D = 10.**tmp[:,1:9] # sample dilution
			k = tmp[:,-8:] # number of positive wells @ each dilution
			dat[f'{expt}_sumVDnmk'] = -pdic['Vinoc']/pdic['S']*numpy.sum(D*(4-k),axis=1)
			dat[f'{expt}_VD'] = -pdic['Vinoc']/pdic['S']*D
			dat[f'{expt}_k'] = k

	else: # LikeSSR
		# --- log10 Vsin in SC
		# Discard points at t=0
		tmp = numpy.loadtxt(data_dir+'/sH1N1_SC_SIN.dat')[3:].T
		# Take ln(10)*log10 V + ln(Vol) to get ln(SIN count)
		lvn = numpy.log(10.)*tmp[1] + numpy.log(pdic['S'])
		# Reshape as 3 columns (one per replicate)
		# 	& then rotate so row=repeats, col=time
		dat['SC_lnVsin'] = numpy.reshape(lvn, (-1,3)).T
		dat['SC_lnVsin_norm'] = -0.5 / (numpy.log(10.) * pdic['sigma_SC'])**2.
		# --- log10 Vsin in MC
		tmp = numpy.loadtxt(data_dir+'/sH1N1_MC_SIN.dat').T
		lvn = numpy.log(10.)*tmp[1] + numpy.log(pdic['S'])
		lvn = numpy.reshape(lvn, (-1,3)).T
		# True where data is nan
		mask = numpy.isnan(lvn)
		# Set nan data = to LLoD & convert units from SIN/ml to ln[SIN]
		lvn[mask] = numpy.log( pdic['LLoD']*pdic['S'] )
		# Array that excludes where SIN was unmeasurable (nan)
		dat['MC_lnVsin'] = numpy.ma.masked_array(lvn,mask)
		# Array that only includes where SIN was unmeasurable (nan)
		dat['MC_lnVsin_nan'] = numpy.ma.masked_array(lvn,mask==False)
		dat['MC_lnVsin_norm'] = -0.5 / (numpy.log(10.) * pdic['sigma_MC'])**2.
	return dat
