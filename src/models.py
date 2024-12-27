#!/usr/bin/env python3
# Copyright (C) 2019-2024
#       Christian Quirouette <https://github.com/cquir>
#       Catherine Beauchemin <cbeau@users.sourceforge.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.	If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import numpy
import scipy

class SMSolution(object):
	""" Solution class for the stochastic model, similar to OdeSolution. """
	def __init__(self):
		# (len(t),) vector of solution times (every tprint hour)
		self.t = []
		# (len(y0),len(t)) array that holds the solution
		self.y = []
		# Set to True if no issue arises
		self.success = False


def solve_stoch_ivp( params, tstart, tend, y, pcrit=0.05, tprint=0.05, meanfield=False ):
	""" Solve the virus infection time-course using the SM model or its mean-field counterpart. """
	# tprint is time interval in hours at which to record data
	# 0.001 means every 3.6 seconds
	# 0.05 means every 3 minutes

	# Pre-compute some stuff
	pdic = params.pardict
	rng = params.rng
	p_over_prna = pdic['rho']/pdic['rho_rna']
	beta_over_S = pdic['beta']/pdic['S']

	# Some pre-calcs for more rapid calculation of time step dt
	#	dt = pcrit / max[ beta/S*T, c, crna, nE/tauE, nI/tauI ]
	# but T changes at each step. We can pre-compute factors to save time
	# dt = pcrit/max(a,b*T) = (pcrit/b)/max(a/b,T) = pcob/max(aob,T)
	pcob = pcrit/beta_over_S
	aob = max(pdic['c'],pdic['crna'],params.spares['kE'],params.spares['dI'])/beta_over_S

	# Continuous populations (deterministic)
	if meanfield:
		rbin = lambda n,p: n*p
		rpoi = lambda lam: lam
		rmul = lambda n,pvals: [n*pvals[0],n*pvals[1],n*(1-pvals[0]-pvals[1])]
	# Discrete populations (stochastic)
	else:
		# Random number generator instance
		rng = params.rng
		rbin = lambda n,p: rng.binomial(n,p)
		rpoi = lambda lam: rng.poisson(lam)
		rmul = lambda n,pvals: rng.multinomial(n,pvals)

	# Extract variables
	(T,E,I,V,Vrna) = (y[0],y[1:1+pdic['nE']],y[-2-pdic['nI']:-2],y[-2],y[-1])
	# Next time at which to print solution for graphing
	tnext = -1.0
	# Array to store graphing solutions
	sol = SMSolution()

	# Starting time
	t = tstart
	# Store value at tstart
	sol.t.append( t )
	sol.y.append( numpy.hstack((T,E,I,V,Vrna)) )

	while t < tend:
		# Determine step size to take based on sensitivity pcrit
		dt = min( pcob/max(aob,T), tend-t )
		if dt == 0.0:
			break
		# Store every "teval" points, and store 1st/last point
		if t >= tnext:
			sol.t.append( t )
			sol.y.append( numpy.hstack((T,E,I,V,Vrna)) )
			tnext = t + tprint
		# Sample how many events occur in time step 
		# Last item in pvals makes sure pvals adds to 1 so value irrelevant
		[Vclr,Vabs,Vleft] = rmul(n=V,pvals=[dt*pdic['c'],dt*beta_over_S*T,0])
		Vnabrt = rbin(n=Vabs,p=pdic['gamma'])
		if (T<1) and meanfield:
			Ninf = 0; T = 0
		elif (T>=1) and meanfield:
			Ninf = T*(1.-((T-1.)/T)**Vnabrt)
		else:
			Ninf = len(numpy.unique(rng.integers(T,size=Vnabrt)))
		EiOut = rbin( n=E, p=dt*params.spares['kE'] )
		IjOut = rbin( n=I, p=dt*params.spares['dI'] )
		Vrnaprod = rpoi( lam=dt*pdic['rho_rna']*numpy.sum(I) )
		Vprod = rbin( n=Vrnaprod, p=p_over_prna )
		Vrnaclr = rbin( n=Vrna, p=dt*pdic['crna'] )
		# Update the system
		t += dt
		Vrna += Vrnaprod - Vabs - Vrnaclr 
		V = Vleft + Vprod
		T -= Ninf
		E[0]  += Ninf - EiOut[0]
		E[1:] -= numpy.diff(EiOut)
		I[0]  += EiOut[-1] - IjOut[0]
		I[1:] -= numpy.diff(IjOut)
	# Store value at tend
	sol.t = numpy.array( sol.t + [t] )
	sol.y = numpy.vstack( sol.y + [numpy.hstack((T,E,I,V,Vrna))] ).T
	# We got to the end, so that's a success ;)
	sol.sucess = True
	return sol


def solve_meanfield_ivp(*args,**kwargs):
	""" Solve the virus infection time-course using the SM's mean-field counterpart. """
	return solve_stoch_ivp(*args, meanfield=True, **kwargs)


def solve_ode_ivp( params, tstart, tend, y ):
	""" Solve the virus infection time-course using the ODE model. """
	def ODE(t,y,params,beta_over_S):
		# Extact variables
		(T,E,I,V,Vrna) = (y[0],y[1:1+pdic['nE']],y[-2-pdic['nI']:-2],y[-2],y[-1])
		# Pre-calculation of shared terms
		sI = numpy.sum(I)
		bTV = beta_over_S*T*V
		kEE = params.spares['kE']*E
		dII = params.spares['dI']*I
		# ODEs
		dT = -pdic['gamma']*bTV
		dE1 = pdic['gamma']*bTV-kEE[0]
		dEi = -numpy.diff(kEE)
		dI1 = kEE[-1]-dII[0]
		dIi = -numpy.diff(dII)
		dV = pdic['rho']*sI - pdic['c']*V - bTV
		dVrna = pdic['rho_rna']*sI - pdic['crna']*Vrna - bTV
		return numpy.hstack((dT,dE1,dEi,dI1,dIi,dV,dVrna))
	pdic = params.pardict
	beta_over_S = pdic['beta']/pdic['S']
	return scipy.integrate.solve_ivp(ODE,[tstart,tend],y,args=(params,beta_over_S),method='BDF')


def infection_simulator(params,tsamples,label,meanfield=False,pcrit=0.05,tprint=0.05):
	""" Reproduces the course of the experimental SC and MC infection as reported in
			https://doi.org/10.1038/srep24154. """
	pdic = params.pardict
	if pdic['SM']:
		# Require running stoch solver in units of IV
		assert pdic['IV']
	# Round initial conditions to nearest int if MM=SM or MFM, not if MM=ODE
	toround = lambda x: x if not pdic['SM'] else round(x)
	# Solver to use (ODE, MFM or SM)
	if not pdic['SM']:
		solve_ivp = solve_ode_ivp
	elif meanfield:
		solve_ivp = solve_meanfield_ivp
	else:
		solve_ivp = solve_stoch_ivp

	# Initialise solution-storing arrays
	tVVrnas = []; V_samples = []; Vrna_samples = []

	# Initial conditions (V=y[-2], Vrna=y[-1])
	y0 = numpy.zeros(pdic['nE']+pdic['nI']+3, dtype=int if (pdic['SM'] and not meanfield) else float)
	# T=y[0] (uninfected cells)
	y0[0] = toround( pdic['Ncells'] )

	# Initializing SC infection
	if label == 'SC':
		# V(-1h)
		y0[-2] = toround( params.spares['Vm1sc'] )
		# Vrna(-1h)
		y0[-1] = toround( pdic['Vm1scrna'] )
		# Don't allow Vrna > V when [V]=IV
		if pdic['IV'] and (y0[-2] > y0[-1]):
			return None
		# Progress SC solution through incubation period (-1h to 0h)
		sol = solve_ivp( params, -1.0, 0.0, y0 )
		# End value becomes initial conditions at t=0 for T,E,I
		y0 = sol.y[:,-1]
		# Set V(0) and Vrna(0) to experimentally measured post-rinse values
		y0[-2] = toround( pdic['V0pr'] / params.spares['PVestED'] )
		y0[-1] = toround( pdic['V0prrna'] )

	# Initializing MC infection
	else: # MC assay
		y0[-2] = toround( pdic['V0mc'] / params.spares['PVestED'] )
		y0[-1] = toround( pdic['V0mcrna'] )
		V_samples.append(y0[-2]); Vrna_samples.append(y0[-1])

	# If V > Vrna when [V]=IV return None because it's unphysical
	if pdic['IV'] and (y0[-2] > y0[-1]):
		return None

	# Now proceed with infection from now
	for i in range(len(tsamples)-1):
		sol = solve_ivp( params, tsamples[i], tsamples[i+1], y0 )
		# Store V and Vrna variables for plotting
		tVVrnas.append( numpy.vstack((sol.t, sol.y[-2:])) )
		# Store variables to be compare against data
		V_samples.append( sol.y[-2][-1] )
		Vrna_samples.append( sol.y[-1][-1] )
		# Initial condition at start of next inter-sampling interval
		y0 = sol.y[:,-1].copy()
		# Apply rinse factor to V and Vrna, and repeat
		if pdic['SM'] and not meanfield:
			# Rinse Vrna
			y0[-1] = params.rng.binomial(n=y0[-1], p=pdic['frinse'])
			# Rinse V but ensure that
			# 	you're not randomly left w/ more V than Vrna
			y0[-2] = min( y0[-1], params.rng.binomial(n=y0[-2], p=pdic['frinse']) )
		else:
			y0[-2] *= pdic['frinse']
			y0[-1] *= pdic['frinse']

	return numpy.hstack(tVVrnas), numpy.array(V_samples), numpy.array(Vrna_samples)
