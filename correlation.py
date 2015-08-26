import numpy as np
import matplotlib.pyplot as plt
import mosquitos_lib as mos
import datos_plot_lib as dpl
from scipy.stats.stats import pearsonr

# time units are in days
dt = 0.01 # time steep (0.01 of a day)
T = 53*7*1.0/dt # how many time steeps as function of number of days (weeks*days)

#defining the numpy needed arrays
m = np.zeros(T)
p = np.zeros(T)
mh = np.zeros(T)
ml = np.zeros(T)
ph = np.zeros(T)
pl = np.zeros(T)
m_qs = np.zeros(T)
mh_qs = np.zeros(T)
ml_qs = np.zeros(T)
cases = np.zeros(T)
#defining numpy arrays for correlation
cases_per_week = np.zeros(int((T)*dt/7))
mosquitos_med = np.zeros(int((T)*dt/7))
mosquitos_high = np.zeros(int((T)*dt/7))
mosquitos_low = np.zeros(int((T)*dt/7))
mosquitos_qs_med = np.zeros(int((T)*dt/7))
mosquitos_qs_high = np.zeros(int((T)*dt/7))
mosquitos_qs_low = np.zeros(int((T)*dt/7))
corr_med = np.zeros(int((T)*dt/7))
corr_high = np.zeros(int((T)*dt/7))
corr_low = np.zeros(int((T)*dt/7))
corr_med_qs = np.zeros(int((T)*dt/7))
corr_high_qs = np.zeros(int((T)*dt/7))
corr_low_qs = np.zeros(int((T)*dt/7))

#defining initial correlation for comparison
corr_med_max = 0.0
corr_high_max = 0.0
corr_low_max = 0.0
corr_med_qs_max = 0.0
corr_high_qs_max = 0.0
corr_low_qs_max = 0.0


#initial conditions
m[0] = 0.5
p[0] = 0.3
mh[0] = 0.5
ph[0] = 0.3
ml[0] = 0.5
pl[0] = 0.3

#initial parameter values
k = 0.5
w = 0.05
epsilon = 0.14*1.05042
phi = 1.0
C = 1.0
pi = 0.12

#initializing the model for med, high and low temperature
mod = mos.model(k,w,epsilon,phi,pi,C)
modh = mos.model(k,w,epsilon,phi,pi,C)
modl = mos.model(k,w,epsilon,phi,pi,C)

#binding temperature data arrays to models
data = dpl.data("DATOS") #loading the data
region = 4 #index of the region
print "region %s" % data.temp_min[region][0]
mod.set_monthly_temperature(data.temp_med[region][1:])
modh.set_monthly_temperature(data.temp_max[region][1:])
modl.set_monthly_temperature(data.temp_min[region][1:])

for i in range(1,int(T)):
	#Set the time steep temperature
	mod.set_model_calendar(i*dt) 
	modh.set_model_calendar(i*dt) 
	modl.set_model_calendar(i*dt)
	#Obtain the quasy stationary solutions
	m_qs[i-1] = mod.m0_qs() # med temperature
	mh_qs[i-1] = modh.m0_qs() # high temperature
	ml_qs[i-1]  = modl.m0_qs() # low temperature
	# Obtain the actual model prediction
	m[i] = m[i-1] + (mod.k*mod.w*p[i-1] - mod.epsilon*m[i-1])*dt
	p[i] = p[i-1] + (mod.phi*m[i-1]*(1.0-p[i-1]/mod.C)-(mod.pi + mod.w)*p[i-1])*dt
	mh[i] = mh[i-1] + (modh.k*modh.w*ph[i-1] - modh.epsilon*mh[i-1])*dt
	ph[i] = ph[i-1] + (modh.phi*mh[i-1]*(1.0-ph[i-1]/modh.C)-(modh.pi + modh.w)*ph[i-1])*dt
	ml[i] = ml[i-1] + (modl.k*modl.w*pl[i-1] - modl.epsilon*ml[i-1])*dt
	pl[i] = pl[i-1] + (modl.phi*ml[i-1]*(1.0-pl[i-1]/modl.C)-(modl.pi + modl.w)*pl[i-1])*dt
	# Obtain the cases from data
	week=int((i)*dt/7)
	week= week - int(week/53)*53
	cases[i] = (1000000.0*data.cases[region][int(week)+1])/data.population[region][1] #percentage of population that got infected that week

#build the point pairs for the correlation
print "Building correlation plots..."
for tau in range(0,53):
	for week in range(1,53+1):
		i=int(((week - 1 - tau)*7 + 3)*(1.0/dt))
		if(i<0):
			i = int(T) + i # left empty part of array due to tau shift is filled with right part (as a circle)
		cases_per_week[week-1] = (1000000.0*data.cases[region][week])/data.population[region][1]
		mosquitos_med[week-1] = m[i]
		mosquitos_high[week-1] = mh[i]
		mosquitos_low[week-1] = ml[i]
		mosquitos_qs_med[week-1] = m_qs[i]
		mosquitos_qs_high[week-1] = mh_qs[i]
		mosquitos_qs_low[week-1] = ml_qs[i]
	#obtain de correlation for each curve with the corresponding shift 'tau'	
	corr_med[tau] = np.corrcoef(mosquitos_med, cases_per_week)[0, 1]
	corr_high[tau] = np.corrcoef(mosquitos_high,cases_per_week)[0,1]
	corr_low[tau] = np.corrcoef(mosquitos_low, cases_per_week)[0,1]
	corr_med_qs[tau] = np.corrcoef(mosquitos_qs_med, cases_per_week)[0,1]
	corr_high_qs[tau] = np.corrcoef(mosquitos_qs_high, cases_per_week)[0,1]
	corr_low_qs[tau] = np.corrcoef(mosquitos_qs_low, cases_per_week)[0,1]
	#find the shift of the max correlation
	if(corr_med_max < corr_med[tau]):
		corr_med_max = corr_med[tau]
		tau_med = tau
	if(corr_high_max < corr_high[tau]):
		corr_high_max = corr_high[tau]
		tau_high = tau
	if(corr_low_max < corr_low[tau]):
		corr_low_max = corr_low[tau]
		tau_low = tau
	if(corr_med_qs_max < corr_med_qs[tau]):
		corr_med_qs_max = corr_med_qs[tau]
		tau_med_qs = tau
	if(corr_high_qs_max < corr_high_qs[tau]):
		corr_high_qs_max = corr_high_qs[tau]
		tau_high_qs = tau
	if(corr_low_qs_max < corr_low_qs[tau]):
		corr_low_qs_max = corr_low_qs[tau]
		tau_low_qs = tau
			
print "corr_med = %f" % corr_med_max
print "tau_med = %f" % tau_med
print "corr_high = %f" % corr_high_max
print "tau_high = %f" % tau_high
print "corr_low = %f" % corr_low_max
print "tau_low = %f" % tau_low
print "corr_med_qs = %f" % corr_med_qs_max
print "tau_med_qs = %f" % tau_med_qs
print "corr_high_qs = %f" % corr_high_qs_max
print "tau_high_qs = %f" % tau_high_qs
print "corr_low_qs = %f" % corr_low_qs_max
print "tau_low_qs = %f" % tau_low_qs 

plt.figure(1)
plt.subplot(2,3,1)
plt.plot(corr_med)
plt.subplot(2,3,2)
plt.plot(corr_high)
plt.subplot(2,3,3)
plt.plot(corr_low)
plt.subplot(2,3,4)
plt.plot(corr_med_qs)
plt.subplot(2,3,5)
plt.plot(corr_high_qs)
plt.subplot(2,3,6)
plt.plot(corr_low_qs)

plt.plot()
plt.show()
#Resultados: Chiapas (4)
#corr_med = 0.562622
#tau_med = 6.000000
#corr_high = 0.759529
#tau_high = 51.000000
#corr_low = 0.714105
#tau_low = 5.000000
#corr_med_qs = 0.577028
#tau_med_qs = 14.000000
#corr_high_qs = 0.772670
#tau_high_qs = 45.000000
#corr_low_qs = 0.681376
#tau_low_qs = 7.000000
