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
			shiftAlerts[i]=(alerts[i] + tau - weeks)
		else:
			if((alerts[i] + tau)<0):
				shiftAlerts[i]=(alerts[i] + tau + weeks)
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
			return count_hits_helper(dtau -1, shift_alerts(leftAlerts,1), outbreaks[negMask],np.append(hits,newHits)) #se comento a partir de la entrada 1008 de la base de datos
			#return count_hits_helper(dtau-1,shift_alerts(alerts,1),outbreaks,np.append(hits,newHits)) #se descomento a partir de 1584 # pruebas a partir de 2160# teoria a partir de 3021

	offsetAlerts = shift_alerts(alerts,tau - dtau/2)
	return count_hits_helper(dtau,offsetAlerts,outbreaks,np.array([]))

def create_region_hits(env):
	def region_hits_function(region):
		alerts = fbt.model_alerts(region,data,delta,c,tauVar)
		outbreaks = fbt.outbreaks(data.cases[region][1:],delta,e,tauVar)
		hits = count_hits(tau,dtau,alerts,outbreaks)[0]
		numHits = len(hits)
		#print "----------------------------"
		#print "Alerts :" 
		#print alerts
		#print "Outbreaks:"
		#print outbreaks
		#print "Hits:"
		#print hits
		#print "-----------------------"
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
	
def insertSQL(db,table,row):
		import sqlite3 as sql
		con = sql.connect(db)
		with con:
			cur = con.cursor()
			cur.execute("INSERT INTO " + table + " VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", row)
	
if __name__ == '__main__':
	env1 = Bunch(data = dpl.data("DATOS"))
	env1.e = Bunch(var = 1.0)
	env1.e.mean = 1.0
	env1.c = Bunch(mean = 1.0)
	env1.c.var = 1.0
	pool = mp.Pool(processes=8)
	# env1.delta = 4 # time lag between two mesoures to compare the increment in mosquitos and incidence  (to stablish treshold)
	#env1.c = 1.96 # number of variances (incidence variance) to stablish warnning treshold
	#env1.tauVar = 4 # period of time to calculate local variance and mean of incidence
	#env1.e = 1.96 # number of variances (incidence variance) to declare outbreak
	#env1.dtau = 3 # interval of time tolerance of possible outbreak predicted by the early warnning
	#env1.tau = 8 # anticipation of the early warnning system
	
	MEC = 1.0 
	FPC = 2.0
	bestEfficiency = 0.0
	lessAsymetricCost = 0.0
	lessSymetricCost = 0.0
	for env1.delta in range(1,2,1):
		#for env1.c.mean in np.arange(0.3,3.6,0.5):
			for env1.c.var in np.arange(0.0,0.1,0.01):
				for env1.tauVar in range(10,11,1):
					for env1.e.var in np.arange(2.5,3.5,0.2):
						for env1.e.mean in np.arange(2.1,3.10,0.2):
							env1.c.mean = env1.e.mean
							for env1.dtau in range(50,51,2):
									for env1.tau in range(0,1):
		
										cumOutbreaks = 0
										cumAlerts = 0
										cumHits = 0
										
										doing = [pool.apply_async(picklable,args=(create_region_hits,env1,region,)) for region in range(32)]
										local_results = np.array([p.get() for p in doing])
										t_local_results = local_results.transpose()
		
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
										
										row = [env1.delta, env1.c.var,env1.tauVar,env1.e.var,env1.dtau,env1.tau,MEC,FPC,efficiency,falseFraction,AsymetricCost,SymetricCost,int(cumHits),int(cumOutbreaks),int(cumAlerts),env1.e.mean,env1.c.mean]
										insertSQL('dengue_results.db','find_best_treshold_mean_mean_diapause',row)
				
										if(efficiency > bestEfficiency ):
											bestEfficiencyRow = list(row)
										if(AsymetricCost < lessAsymetricCost):
											lessAsymetricCostRow = list(row)
										if(SymetricCost < lessSymetricCost):
											lessSymetricCostRow = list(row)
			
	#insertSQL('dengue_results.db','best_treshold',bestEfficiencyRow)
	#insertSQL('dengue_results.db','best_treshold',lessAsymetricCostRow)
	#insertSQL('dengue_results.db','best_treshold',lessSymetricCostRow)
