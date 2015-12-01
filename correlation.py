import numpy as np
import matplotlib.pyplot as plt
import mosquitos_lib as mos
import datos_plot_lib as dpl
from scipy.stats.stats import pearsonr
import scipy, scipy.stats

import matplotlib

matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18)
matplotlib.rc('axes', labelsize=18)

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

corr_med_temp = np.zeros(int((T)*dt/7))
corr_high_temp = np.zeros(int((T)*dt/7))
corr_low_temp = np.zeros(int((T)*dt/7))
shifted_temperature_min = np.zeros(int((T)*dt/7))
shifted_temperature_max = np.zeros(int((T)*dt/7))
shifted_temperature_med = np.zeros(int((T)*dt/7))

#defining initial correlation for comparison
corr_med_max = 0.0
corr_high_max = 0.0
corr_low_max = 0.0
corr_med_qs_max = 0.0
corr_high_qs_max = 0.0
corr_low_qs_max = 0.0
corr_med_temp_max = 0.0
corr_high_temp_max = 0.0
corr_low_temp_max = 0.0

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
region = 2 #25 #index of the region
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
	#if(week < 39):
	#	ml[i]=0.01
	#	pl[i]=0.01
	
#####

temp = []
for i in range(0,len(data.cases[region])-1):
	temp.append(data.weekly_temperature(i,region,"min"))

ShiftTemp = np.zeros(len(temp))
ShiftModel = np.zeros(len(temp))


for i in range(0,len(temp)):
    weekR = i - 9
    if(weekR < 0):
        weekR = len(temp) + weekR
    ShiftTemp[i] = temp[weekR]
    
fig2, ax1 = plt.subplots()

ax1.plot(ShiftTemp,data.cases[region][1:],'bo')
ax1.set_xlabel('Temperature (C)',)
ax1.set_ylabel('Casos')
#plt.show()

#ajuste a una cuadratica
import scipy.optimize as optimization

def func(x, a, b, c):
    return a + b*x + c*x*x
    
def func3(x, a, b, c, d):
    return a + b*x + c*x*x + d*x*x*x

xdata = ShiftTemp
ydata = data.cases[region][1:]
x0 = np.array([10.0, 20.0, 1.0])
sigma = np.ones(len(xdata))

print optimization.curve_fit(func, xdata, ydata, x0, sigma)

fig3, ax1 = plt.subplots()
ax1.plot(xdata,func(xdata,741.195,-99.679,3.277),'g')
ax1.plot(ShiftTemp,data.cases[region][1:],'bo')
ax1.set_xlabel('Temperature (C)',)
ax1.set_ylabel('Casos')
ax1.axis([12, 25, 0, 400])
#plt.show()

#####


#### Model prediction
prediction = np.zeros(len(data.cases[region])-1)
modelPrediction = np.zeros(len(data.cases[region])-1)
modelPrediction3 = np.zeros(len(data.cases[region])-1)
for i in range(0,len(data.cases[region])-1):
    prediction[i]=func3(ShiftTemp[i], -1.90294525e+03,   3.72969314e+02,  -2.38602639e+01,
         5.01366047e-01)
    modelPrediction[i] = func(mosquitos_low[i],17.08754532, -147.0048907 ,  178.19032476)
    modelPrediction3[i] = func3(mosquitos_low[i], 0.28796529,  121.13014887, -250.92946308,  164.49342213)

fig5, ax1 = plt.subplots()
ax1.plot(range(0,len(data.cases[region])-1),data.cases[region][1:],'go')
ax1.plot(range(0,len(data.cases[region])-1),prediction,'r')
#ax1.plot(range(0,len(data.cases[region])-1),modelPrediction,'b')
ax1.plot(range(0,len(data.cases[region])-1),modelPrediction3,'b')
ax1.axis([0, 53, 0, 450])
ax1.set_xlabel('Time (weeks)',)
ax1.set_ylabel('Incidence')
fig5.savefig('MixedModel3.pdf', bbox_inches='tight')
plt.show()

ChiT=0.0
ChiM2 = 0.0
ChiM3 = 0.0
for i in range(1,len(data.cases[region][1:])):
		ChiT = ChiT + (data.cases[region][i] - prediction[i])*(data.cases[region][i] - prediction[i])
		ChiM2 = ChiM2 + (data.cases[region][i] - modelPrediction[i])*(data.cases[region][i] - modelPrediction[i])
		ChiM3 = ChiM3 + (data.cases[region][i] - modelPrediction3[i])*(data.cases[region][i] - modelPrediction3[i])
		
print 'ChiT= %f' % ChiT
print 'ChiM2= %f' % ChiM2
print 'ChiM3= %f' % ChiM3

import scipy, scipy.stats
 
observed_values=scipy.array(data.cases[region][1:])
expected_values=scipy.array(prediction)
 
print scipy.stats.chisquare(observed_values, f_exp=expected_values)

expected_values=scipy.array(modelPrediction)
print scipy.stats.chisquare(observed_values, f_exp=expected_values)

expected_values=scipy.array(modelPrediction3)
print scipy.stats.chisquare(observed_values, f_exp=expected_values)
#####################################################################	
#2758.942701313
#3485.203593469
#3218.379379688
#######
	
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
		shifted_temperature_min[week-1] = data.weekly_temperature(week - tau, region , "min")
		shifted_temperature_max[week-1] = data.weekly_temperature(week - tau, region , "max")
		shifted_temperature_med[week-1] = data.weekly_temperature(week - tau, region , "med")
	#obtain de correlation for each curve with the corresponding shift 'tau'	
	corr_med[tau] = np.corrcoef(mosquitos_med, cases_per_week)[0, 1]
	corr_high[tau] = np.corrcoef(mosquitos_high,cases_per_week)[0,1]
	corr_low[tau] = np.corrcoef(mosquitos_low, cases_per_week)[0,1]
	corr_med_qs[tau] = np.corrcoef(mosquitos_qs_med, cases_per_week)[0,1]
	corr_high_qs[tau] = np.corrcoef(mosquitos_qs_high, cases_per_week)[0,1]
	corr_low_qs[tau] = np.corrcoef(mosquitos_qs_low, cases_per_week)[0,1]
	
	corr_med_temp[tau] = np.corrcoef(shifted_temperature_med, cases_per_week)[0,1]
	corr_high_temp[tau] = np.corrcoef(shifted_temperature_max, cases_per_week)[0,1]
	corr_low_temp[tau] = np.corrcoef(shifted_temperature_min, cases_per_week)[0,1]
	
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
	if(corr_med_temp_max < corr_med_temp[tau]):
		corr_med_temp_max = corr_med_temp[tau]
		tau_med_temp = tau
	if(corr_high_temp_max < corr_high_temp[tau]):
		corr_high_temp_max = corr_high_temp[tau]
		tau_high_temp = tau
	if(corr_low_temp_max < corr_low_temp[tau]):
		corr_low_temp_max = corr_low_temp[tau]
		tau_low_temp = tau
			
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

print "corr_med_temp = %f" % corr_med_temp_max
print "tau_med_temp = %f" % tau_med_temp
print "corr_high_temp = %f" % corr_high_temp_max
print "tau_high_temp = %f" % tau_high_temp
print "corr_low_temp = %f" % corr_low_temp_max
print "tau_low_temp = %f" % tau_low_temp 

plt.figure(1)
plt.subplot(2,3,1)
plt.plot(corr_med)
plt.plot(corr_med_temp,'--')
plt.subplot(2,3,2)
plt.plot(corr_high)
plt.plot(corr_high_temp,'--')
plt.subplot(2,3,3)
plt.plot(corr_low)
plt.plot(corr_low_temp,'--')
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
#corr_med_temp = 0.564111
#tau_med_temp = 13.000000
#corr_high_temp = 0.522739
#tau_high_temp = 18.000000
#corr_low_temp = 0.681588 ---
#tau_low_temp = 8.000000

#region Baja California Sur (2)
#corr_med = 0.595154
#tau_med = 4.000000
#corr_high = 0.730854
#tau_high = 29.000000
#corr_low = 0.835484
#tau_low = 4.000000
#corr_med_qs = 0.559890
#tau_med_qs = 15.000000
#corr_high_qs = 0.677147
#tau_high_qs = 36.000000
#corr_low_qs = 0.184277
#tau_low_qs = 41.000000
#corr_med_temp = 0.764122
#tau_med_temp = 10.000000
#corr_high_temp = 0.730422
#tau_high_temp = 11.000000
#corr_low_temp = 0.786509 ---
#tau_low_temp = 9.000000

#region Guerrero (12)----
#corr_med = 0.481825
#tau_med = 10.000000
#corr_high = 0.543532
#tau_high = 50.000000
#corr_low = 0.589716
#tau_low = 5.000000
#corr_med_qs = 0.589457
#tau_med_qs = 13.000000
#corr_high_qs = 0.584971
#tau_high_qs = 46.000000
#corr_low_qs = 0.584231
#tau_low_qs = 9.000000
#corr_med_temp = 0.580331
#tau_med_temp = 13.000000
#corr_high_temp = 0.576175
#tau_high_temp = 17.000000
#corr_low_temp = 0.579712 ---
#tau_low_temp = 10.000000

#region Sonora (25)
#corr_med = 0.473598
#tau_med = 8.000000
#corr_high = 0.533123
#tau_high = 35.000000
#corr_low = 0.738731
#tau_low = 14.000000
#corr_med_qs = 0.469121
#tau_med_qs = 11.000000
#corr_high_qs = 0.454511
#tau_high_qs = 39.000000
#corr_low_qs = 0.370419
#tau_low_qs = 31.000000
#corr_med_temp = 0.495189
#tau_med_temp = 20.000000
#corr_high_temp = 0.530424 ---
#tau_high_temp = 22.000000
#corr_low_temp = 0.446617 
#tau_low_temp = 19.000000

#region Sonora (just 2 peak)
#corr_med = 0.473598
#tau_med = 8.000000
#corr_high = 0.533123
#tau_high = 35.000000
#corr_low = 0.841084
#tau_low = 2.000000
#corr_med_qs = 0.469121
#tau_med_qs = 11.000000
#corr_high_qs = 0.454511
#tau_high_qs = 39.000000
#corr_low_qs = 0.370419
#tau_low_qs = 31.000000
#corr_med_temp = 0.495189
#tau_med_temp = 20.000000
#corr_high_temp = 0.530424
#tau_high_temp = 22.000000
#corr_low_temp = 0.446617
#tau_low_temp = 19.000000

#region Colima (7)
#corr_med = 0.491214
#tau_med = 52.000000
#corr_high = 0.610621
#tau_high = 44.000000
#corr_low = 0.689966
#tau_low = 5.000000
#corr_med_qs = 0.642413
#tau_med_qs = 7.000000
#corr_high_qs = 0.702156
#tau_high_qs = 45.000000
#corr_low_qs = 0.728621
#tau_low_qs = 6.000000
#corr_med_temp = 0.676070
#tau_med_temp = 9.000000
#corr_high_temp = 0.642097
#tau_high_temp = 18.000000
#corr_low_temp = 0.724862 ---
#tau_low_temp = 7.000000

#region Oaxaca (19)
#corr_med = 0.563627
#tau_med = 9.000000
#corr_high = 0.702489
#tau_high = 46.000000
#corr_low = 0.678785
#tau_low = 4.000000
#corr_med_qs = 0.581415
#tau_med_qs = 15.000000
#corr_high_qs = 0.769232
#tau_high_qs = 42.000000
#corr_low_qs = 0.603069
#tau_low_qs = 9.000000
#corr_med_temp = 0.564582
#tau_med_temp = 16.000000
#corr_high_temp = 0.619821
#tau_high_temp = 21.000000
#corr_low_temp = 0.606476 ---
#tau_low_temp = 10.000000

#region Sinaloa (24)
#corr_med = 0.557325
#tau_med = 1.000000
#corr_high = 0.640091
#tau_high = 39.000000
#corr_low = 0.646854
#tau_low = 5.000000
#corr_med_qs = 0.580290
#tau_med_qs = 5.000000
#corr_high_qs = 0.709671
#tau_high_qs = 41.000000
#corr_low_qs = 0.291471
#tau_low_qs = 51.000000
#corr_med_temp = 0.507262
#tau_med_temp = 12.000000
#corr_high_temp = 0.584191 ---
#tau_high_temp = 20.000000
#corr_low_temp = 0.573519
#tau_low_temp = 9.000000

