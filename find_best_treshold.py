import numpy as np
import matplotlib.pyplot as plt
import datos_plot_lib as dpl
from scipy.stats.stats import pearsonr
import scipy, scipy.stats
import find_best_treshold_lib as fbt
import matplotlib
import multiprocessing as mp

class Bunch:
	def __init__(self, **kwds):
		self.__dict__.update(kwds)

matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18)
matplotlib.rc('axes', labelsize=18)


def shift_alerts(alerts,tau):
	weeks = 53
	shiftAlerts = np.zeros(len(alerts))
	for i in range(len(alerts)):
		if((alerts[i]+tau)>=weeks):
			shiftAlerts[i]=alerts[i]+tau - weeks
		else:
			shiftAlerts[i]=alerts[i]+tau
	return shiftAlerts

def count_hits(tau,dtau,alerts,outbreaks):	
	def count_hits_helper(dtau, alerts, outbreaks,hits):
		if(dtau < 0):
			return np.array([hits,outbreaks,alerts])
		else:
			mask = np.in1d(outbreaks,alerts)
			newHits = outbreaks[mask]
			negMask = np.logical_not(mask)
			leftAlerts = alerts[np.logical_not(np.in1d(alerts,newHits))]
			return count_hits_helper(dtau -1, shift_alerts(leftAlerts,1), outbreaks[negMask],np.append(hits,newHits))

	offsetAlerts = shift_alerts(alerts,tau)
	return count_hits_helper(dtau,offsetAlerts,outbreaks,np.array([]))

def create_region_hits(env):
	def region_hits_function(region):
		alerts = fbt.model_alerts(region,data,delta,c,tauVar)
		outbreaks = fbt.outbreaks(data.cases[region][1:],delta,e,tauVar)
		numHits = len(count_hits(tau,dtau,alerts,outbreaks)[0])
		return np.array([numHits,len(outbreaks),len(alerts)])
	
	#time units are in weeks
	data = env.data
	delta = env.delta # time lag between two mesoures to compare the increment in mosquitos and incidence  (to stablish treshold)
	c = env.c # number of variances (incidence variance) to stablish warnning treshold
	tauVar = env.tauVar # period of time to calculate local variance and mean of incidence
	e = env.e # number of variances (incidence variance) to declare outbreak

	###
	dtau = env.dtau # interval of time of possible outbreak predicted by the early warnning
	tau = env.tau # anticipation of the early warnning system
	###
	
	return region_hits_function
	
def picklable(creation,env,region):
	return creation(env)(region)
	
if __name__ == '__main__':
	env1 = Bunch(data = dpl.data("DATOS"))
	env1.delta = 4 # time lag between two mesoures to compare the increment in mosquitos and incidence  (to stablish treshold)
	env1.c = 1.96 # number of variances (incidence variance) to stablish warnning treshold
	env1.tauVar = 4 # period of time to calculate local variance and mean of incidence
	env1.e = 1.96 # number of variances (incidence variance) to declare outbreak

		###
	env1.dtau = 3 # interval of time of possible outbreak predicted by the early warnning
	env1.tau = 8 # anticipation of the early warnning system

	MEC = 1.0 
	FPC = 2.0

	cumOutbreaks = 0
	cumAlerts = 0
	cumHits = 0

	pool = mp.Pool(processes=4)

	#region_hits = create_region_hits(env1)
	

	doing = [pool.apply_async(picklable,args=(create_region_hits,env1,region,)) for region in range(32)]
	#local_results = np.array(pool.map(picklable, range(32)))
	#local_results = np.array(map(region_hits,range(32)))
	local_results = np.array([p.get() for p in doing])
	t_local_results = local_results.transpose()
	print local_results

	results = [np.sum(t_local_results[n]) for n in range(len(local_results[0]))] 
	
	cumOutbreaks = cumOutbreaks + results[1]
	cumAlerts = cumAlerts + results[2]
	cumHits = cumHits + results[0]

	if(cumOutbreaks != 0):
		efficiency = float(cumHits)/float(cumOutbreaks)
	else:
		efficiency = -1.0
		print "no good outbreak treshold"

	if(cumAlerts != 0):
		falseFraction = float(cumAlerts - cumHits)/float(cumAlerts)
	else:
		falseFraction = -1.0
		print "no good trigger treshold"
		
	AsymetricCost = MEC*(cumOutbreaks - cumHits) + FPC*(cumAlerts - cumHits)
	SymetricCost = FPC*(cumOutbreaks + cumAlerts - 2.0*cumHits)

	print efficiency
	print falseFraction
	print AsymetricCost
	print SymetricCost
	print cumHits
	print cumOutbreaks
	print cumAlerts
