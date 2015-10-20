import numpy as np
import matplotlib.pyplot as plt
import mosquitos_lib as mos
import datos_plot_lib as dpl
import matplotlib

matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18)
matplotlib.rc('axes', labelsize=18)
#plt.rcParams.update({'axes.labelsize': 'large'})
#50*50
dt = 0.01
T = 430.0*1.0/dt

m = np.zeros(T)
p = np.zeros(T)
mh = np.zeros(T)
ml = np.zeros(T)
ph = np.zeros(T)
pl = np.zeros(T)
m0_qs = np.zeros(T)
p0_qs = np.zeros(T)
m1_qs = np.zeros(T)
cases = np.zeros(T)
temp = np.zeros(T+1)
dataTemp = np.zeros(T+1)
cum_cases = np.zeros(T)
m_cum = np.zeros(T)
calendar = np.zeros(T)
xaxis = np.zeros(T)

#for region in range(0,31):
m[0] = 0.5
p[0] = 0.3
mh[0] = 0.5
ph[0] = 0.3
ml[0] = 0.5
pl[0] = 0.3

cum_cases[0] = 0.0
m_cum[0]= 0.0

k = 0.5
w = 0.05
epsilon = 0.14*1.05042
phi = 1.0
C = 1.0
pi = 0.12

#Fig 1
region = 30 #2
mod = mos.model(k,w,epsilon,phi,pi,C)
modh = mos.model(k,w,epsilon,phi,pi,C)
modl = mos.model(k,w,epsilon,phi,pi,C)
data = dpl.data("DATOS")
mod.set_monthly_temperature(data.temp_med[region][1:])
modh.set_monthly_temperature(data.temp_max[region][1:])
modl.set_monthly_temperature(data.temp_min[region][1:])
print "region %s" % data.temp_min[region][0]
print "poblacion de %s" % data.population[region][0]


for i in range(1,int(T)):
	mod.set_model_calendar(i*dt)
	modh.set_model_calendar(i*dt)
	modl.set_model_calendar(i*dt)
	m0_qs[i-1] = mod.m0_qs()
	m1_qs[i-1] = -mod.Dm0_qs((i-1)*dt,28)/mod.epsilon
	m[i] = m[i-1] + (mod.k*mod.w*p[i-1] - mod.epsilon*m[i-1])*dt
	p[i] = p[i-1] + (mod.phi*m[i-1]*(1.0-p[i-1]/mod.C)-(mod.pi + mod.w)*p[i-1])*dt
	mh[i] = mh[i-1] + (modh.k*modh.w*ph[i-1] - modh.epsilon*mh[i-1])*dt
	ph[i] = ph[i-1] + (modh.phi*mh[i-1]*(1.0-ph[i-1]/modh.C)-(modh.pi + modh.w)*ph[i-1])*dt
	ml[i] = ml[i-1] + (modl.k*modl.w*pl[i-1] - modl.epsilon*ml[i-1])*dt
	pl[i] = pl[i-1] + (modl.phi*ml[i-1]*(1.0-pl[i-1]/modl.C)-(modl.pi + modl.w)*pl[i-1])*dt
	if(mh[i] <= 0.01 and ph[i] <= 0.01):
		mh[i]=0.01
		
	week=int((i)*dt/7)
	calendar[i]=week
	xaxis[i]=float(i)*dt/7.0
	week= week - int(week/53)*53
	cases[i] = data.cases[region][int(week)+1]*(1000.0/data.population[region][1])
	cum_cases[i] = cum_cases[i-1] + cases[i]*dt
	m_cum[i] = m_cum[i-1] + ml[i]*dt
	if(week <= 40):
		mh[i]=0.01
		ph[i]=0.01

###Ploting
fig1, ax1 = plt.subplots()

ax1.plot(xaxis,m,'y--',xaxis,ml,'b--',xaxis,mh,'r--')
ax1.set_xlabel('time (weeks)',)
ax1.set_ylabel('mosquitoes/C')
ax1.axis([0, 53, 0, 2.5])

ax2 = ax1.twinx()
ax2.plot(calendar,cases,'go')

ax2.set_ylabel('cases per 1000 people')

plt.xlim([0, 53])
#fig1.suptitle(data.population[region][0], fontsize=20)
fig1.show()
#fig1.savefig('fig6.pdf', bbox_inches='tight')
## Plot Temperature
fig1, ax1 = plt.subplots()

ax1.set_xlabel('time (weeks)')
ax1.set_ylabel('temperature (C)')
ax1.plot(range(0,53),np.array(data.weekly_temperature_array(region,'min')),'b')
ax1.plot(range(0,53),np.array(data.weekly_temperature_array(region,'max')),'r')
ax1.plot(range(0,53),np.array(data.weekly_temperature_array(region,'med')),'y')

ax2 = ax1.twinx()
ax2.plot(calendar,cases,'go')

ax2.set_ylabel('cases per 1000 people')

plt.xlim([0, 53])
#fig1.suptitle(data.population[region][0], fontsize=20)
fig1.show()
#fig1.savefig('fig6temp.pdf', bbox_inches='tight')
#plt.subplot(2,2,1)
#plt.plot(m,'y--', cases,'g--',ml,'b--',mh,'r--')
#plt.plot(m_cum,'r--', cum_cases/200.0,'g--')
#plt.axis([100, T, 0, 3])	
#FIg 2 ##############################################################
region = 31 #12 #24 #7 #4
mod = mos.model(k,w,epsilon,phi,pi,C)
modh = mos.model(k,w,epsilon,phi,pi,C)
modl = mos.model(k,w,epsilon,phi,pi,C)
data = dpl.data("DATOS")
mod.set_monthly_temperature(data.temp_med[region][1:])
modh.set_monthly_temperature(data.temp_max[region][1:])
modl.set_monthly_temperature(data.temp_min[region][1:])
print "region %s" % data.temp_min[region][0] 
print "poblacion de %s" % data.population[region][0]

for i in range(1,int(T)):
	mod.set_model_calendar(i*dt)
	modh.set_model_calendar(i*dt)
	modl.set_model_calendar(i*dt)
	m0_qs[i-1] = mod.m0_qs()
	m1_qs[i-1] = -mod.Dm0_qs((i-1)*dt,28)/mod.epsilon
	m[i] = m[i-1] + (mod.k*mod.w*p[i-1] - mod.epsilon*m[i-1])*dt
	p[i] = p[i-1] + (mod.phi*m[i-1]*(1.0-p[i-1]/mod.C)-(mod.pi + mod.w)*p[i-1])*dt
	mh[i] = mh[i-1] + (modh.k*modh.w*ph[i-1] - modh.epsilon*mh[i-1])*dt
	ph[i] = ph[i-1] + (modh.phi*mh[i-1]*(1.0-ph[i-1]/modh.C)-(modh.pi + modh.w)*ph[i-1])*dt
	ml[i] = ml[i-1] + (modl.k*modl.w*pl[i-1] - modl.epsilon*ml[i-1])*dt
	pl[i] = pl[i-1] + (modl.phi*ml[i-1]*(1.0-pl[i-1]/modl.C)-(modl.pi + modl.w)*pl[i-1])*dt
	week=int((i)*dt/7)
	week= week - int(week/53)*53
	calendar[i]=week
	xaxis[i]=float(i)*dt/7.0
	cases[i] = data.cases[region][int(week)+1]*(10000.0/data.population[region][1])
	cum_cases[i] = cum_cases[i-1] + cases[i]*dt
	m_cum[i] = m_cum[i-1] + ml[i]*dt

fig2, ax1 = plt.subplots()

ax1.plot(xaxis,m,'y--',xaxis,ml,'b--',xaxis,mh,'r--')
ax1.set_xlabel('time (weeks)')
ax1.set_ylabel('mosquitoes/C')
ax1.axis([0, 53, 0, 2.5])

ax2 = ax1.twinx()
ax2.plot(calendar,cases,'go')
ax2.set_ylabel('cases per 10000 people')

plt.xlim([0, 53])
#fig2.suptitle(data.population[region][0], fontsize=20)
fig2.show()
#fig2.savefig('fig2.pdf', bbox_inches='tight')
## Plot Temperature
fig2, ax1 = plt.subplots()

ax1.set_xlabel('time (weeks)')
ax1.set_ylabel('temperature (C)')
ax1.plot(range(0,53),np.array(data.weekly_temperature_array(region,'min')),'b')
ax1.plot(range(0,53),np.array(data.weekly_temperature_array(region,'max')),'r')
ax1.plot(range(0,53),np.array(data.weekly_temperature_array(region,'med')),'y')

ax2 = ax1.twinx()
ax2.plot(calendar,cases,'go')

ax2.set_ylabel('cases per 10000 people')

plt.xlim([0, 53])
#fig2.suptitle(data.population[region][0], fontsize=20)
#fig2.show()
fig2.savefig('fig2temp.pdf', bbox_inches='tight')

#exit()
#plt.subplot(2,2,2)
#plt.plot(m,'y--', cases,'g--',ml,'b--',mh,'r--')
#plt.plot(cum_cases,'g--',m_cum,'b--')
input()
################################################################
#Fig 3 ##############################################################
region = 32 #24 #25 #5
mod = mos.model(k,w,epsilon,phi,pi,C)
modh = mos.model(k,w,epsilon,phi,pi,C)
modl = mos.model(k,w,epsilon,phi,pi,C)
data = dpl.data("DATOS")
mod.set_monthly_temperature(data.temp_med[region][1:])
modh.set_monthly_temperature(data.temp_max[region][1:])
modl.set_monthly_temperature(data.temp_min[region][1:])
print "region %s" % data.temp_min[region][0] 
print "poblacion de %s" % data.population[region][0]

for i in range(1,int(T)):
	mod.set_model_calendar(i*dt)
	modh.set_model_calendar(i*dt)
	modl.set_model_calendar(i*dt)
	m0_qs[i-1] = mod.m0_qs()
	m1_qs[i-1] = -mod.Dm0_qs((i-1)*dt,28)/mod.epsilon
	m[i] = m[i-1] + (mod.k*mod.w*p[i-1] - mod.epsilon*m[i-1])*dt
	p[i] = p[i-1] + (mod.phi*m[i-1]*(1.0-p[i-1]/mod.C)-(mod.pi + mod.w)*p[i-1])*dt
	mh[i] = mh[i-1] + (modh.k*modh.w*ph[i-1] - modh.epsilon*mh[i-1])*dt
	ph[i] = ph[i-1] + (modh.phi*mh[i-1]*(1.0-ph[i-1]/modh.C)-(modh.pi + modh.w)*ph[i-1])*dt
	ml[i] = ml[i-1] + (modl.k*modl.w*pl[i-1] - modl.epsilon*ml[i-1])*dt
	pl[i] = pl[i-1] + (modl.phi*ml[i-1]*(1.0-pl[i-1]/modl.C)-(modl.pi + modl.w)*pl[i-1])*dt
	week=int((i)*dt/7)
	week= week - int(week/53)*53
	calendar[i]=week
	xaxis[i]=float(i)*dt/7.0
	cases[i] = float(data.cases[region][int(week)+1])*(10000.0/float(data.population[region][1]))
	cum_cases[i] = cum_cases[i-1] + cases[i]*dt
	m_cum[i] = m_cum[i-1] + ml[i]*dt


fig3, ax1 = plt.subplots()

ax1.plot(xaxis,m,'y--',xaxis,ml,'b--',xaxis,mh,'r--')
ax1.set_xlabel('time (weeks)')
ax1.set_ylabel('mosquitoes/C')
ax1.axis([0, 53, 0, 2.5])

ax2 = ax1.twinx()
ax2.plot(calendar,cases,'go')
ax2.set_ylabel('cases per 10000 people')

plt.xlim([0, 53])
#fig3.suptitle(data.population[region][0], fontsize=20)
fig3.show()
#fig3.savefig('fig3.pdf', bbox_inches='tight')
## Plot Temperature
fig3, ax1 = plt.subplots()

ax1.set_xlabel('time (weeks)')
ax1.set_ylabel('temperature (C)')
ax1.plot(range(0,53),np.array(data.weekly_temperature_array(region,'min')),'b')
ax1.plot(range(0,53),np.array(data.weekly_temperature_array(region,'max')),'r')
ax1.plot(range(0,53),np.array(data.weekly_temperature_array(region,'med')),'y')

ax2 = ax1.twinx()
ax2.plot(calendar,cases,'go')

ax2.set_ylabel('cases per 10000 people')

plt.xlim([0, 53])
#fig3.suptitle(data.population[region][0], fontsize=20)
fig3.show()
#fig3.savefig('fig3temp.pdf', bbox_inches='tight')

#plt.subplot(2,2,3)
#plt.plot(m,'y--', cases,'g--',ml,'b--',mh,'r--')
#plt.plot(cum_cases,'g--',m_cum,'b--')
################################################################
#Fig 4##############################################################
region = 33 #7 #29 #6
mod = mos.model(k,w,epsilon,phi,pi,C)
modh = mos.model(k,w,epsilon,phi,pi,C)
modl = mos.model(k,w,epsilon,phi,pi,C)
data = dpl.data("DATOS")
mod.set_monthly_temperature(data.temp_med[region][1:])
modh.set_monthly_temperature(data.temp_max[region][1:])
modl.set_monthly_temperature(data.temp_min[region][1:])
print "region %s" % data.temp_min[region][0]
print "poblacion de %s" % data.population[region][0]


for i in range(1,int(T)):
	mod.set_model_calendar(i*dt)
	modh.set_model_calendar(i*dt)
	modl.set_model_calendar(i*dt)
	m0_qs[i-1] = mod.m0_qs()
	m1_qs[i-1] = -mod.Dm0_qs((i-1)*dt,28)/mod.epsilon
	m[i] = m[i-1] + (mod.k*mod.w*p[i-1] - mod.epsilon*m[i-1])*dt
	p[i] = p[i-1] + (mod.phi*m[i-1]*(1.0-p[i-1]/mod.C)-(mod.pi + mod.w)*p[i-1])*dt
	mh[i] = mh[i-1] + (modh.k*modh.w*ph[i-1] - modh.epsilon*mh[i-1])*dt
	ph[i] = ph[i-1] + (modh.phi*mh[i-1]*(1.0-ph[i-1]/modh.C)-(modh.pi + modh.w)*ph[i-1])*dt
	ml[i] = ml[i-1] + (modl.k*modl.w*pl[i-1] - modl.epsilon*ml[i-1])*dt
	pl[i] = pl[i-1] + (modl.phi*ml[i-1]*(1.0-pl[i-1]/modl.C)-(modl.pi + modl.w)*pl[i-1])*dt
	week=int((i)*dt/7)
	week= week - int(week/53)*53
	calendar[i]=week
	xaxis[i]=float(i)*dt/7.0
	cases[i] = data.cases[region][int(week)+1]*(10000.0/data.population[region][1])
	cum_cases[i] = cum_cases[i-1] + cases[i]*dt
	m_cum[i] = m_cum[i-1] + ml[i]*dt

fig4, ax1 = plt.subplots()

ax1.plot(xaxis,m,'y--',xaxis,ml,'b--',xaxis,mh,'r--')
ax1.set_xlabel('time (weeks)')
ax1.set_ylabel('mosquitoes/C')
ax1.axis([0, 53, 0, 2.5])

ax2 = ax1.twinx()
ax2.plot(calendar,cases,'go')
ax2.set_ylabel('cases per 10000 people')

plt.xlim([0, 53])
#fig4.suptitle(data.population[region][0], fontsize=20)
fig4.show()
#fig4.savefig('fig4.pdf', bbox_inches='tight')
## Plot Temperature
fig4, ax1 = plt.subplots()

ax1.set_xlabel('time (weeks)')
ax1.set_ylabel('temperature (C)')
ax1.plot(range(0,53),np.array(data.weekly_temperature_array(region,'min')),'b')
ax1.plot(range(0,53),np.array(data.weekly_temperature_array(region,'max')),'r')
ax1.plot(range(0,53),np.array(data.weekly_temperature_array(region,'med')),'y')

ax2 = ax1.twinx()
ax2.plot(calendar,cases,'go')

ax2.set_ylabel('cases per 10000 people')

plt.xlim([0, 53])
#fig4.suptitle(data.population[region][0], fontsize=20)
fig4.show()
#fig4.savefig('fig4temp.pdf', bbox_inches='tight')

#plt.subplot(2,2,4)
#plt.plot(m,'y--', cases,'g--',ml,'b--',mh,'r--')
#plt.plot(cum_cases,'g--',m_cum,'b--')
################################################################
#FIg 5 ##############################################################
region = 24 #7 #4
mod = mos.model(k,w,epsilon,phi,pi,C)
modh = mos.model(k,w,epsilon,phi,pi,C)
modl = mos.model(k,w,epsilon,phi,pi,C)
data = dpl.data("DATOS")
mod.set_monthly_temperature(data.temp_med[region][1:])
modh.set_monthly_temperature(data.temp_max[region][1:])
modl.set_monthly_temperature(data.temp_min[region][1:])
print "region %s" % data.temp_min[region][0] 
print "poblacion de %s" % data.population[region][0]

for i in range(1,int(T)):
	mod.set_model_calendar(i*dt)
	modh.set_model_calendar(i*dt)
	modl.set_model_calendar(i*dt)
	m0_qs[i-1] = mod.m0_qs()
	m1_qs[i-1] = -mod.Dm0_qs((i-1)*dt,28)/mod.epsilon
	m[i] = m[i-1] + (mod.k*mod.w*p[i-1] - mod.epsilon*m[i-1])*dt
	p[i] = p[i-1] + (mod.phi*m[i-1]*(1.0-p[i-1]/mod.C)-(mod.pi + mod.w)*p[i-1])*dt
	mh[i] = mh[i-1] + (modh.k*modh.w*ph[i-1] - modh.epsilon*mh[i-1])*dt
	ph[i] = ph[i-1] + (modh.phi*mh[i-1]*(1.0-ph[i-1]/modh.C)-(modh.pi + modh.w)*ph[i-1])*dt
	ml[i] = ml[i-1] + (modl.k*modl.w*pl[i-1] - modl.epsilon*ml[i-1])*dt
	pl[i] = pl[i-1] + (modl.phi*ml[i-1]*(1.0-pl[i-1]/modl.C)-(modl.pi + modl.w)*pl[i-1])*dt
	week=int((i)*dt/7)
	week= week - int(week/53)*53
	calendar[i]=week
	xaxis[i]=float(i)*dt/7.0
	cases[i] = data.cases[region][int(week)+1]*(10000.0/data.population[region][1])
	cum_cases[i] = cum_cases[i-1] + cases[i]*dt
	m_cum[i] = m_cum[i-1] + ml[i]*dt

fig5, ax1 = plt.subplots()

ax1.plot(xaxis,m,'y--',xaxis,ml,'b--',xaxis,mh,'r--')
ax1.set_xlabel('time (weeks)')
ax1.set_ylabel('density/C')
ax1.axis([0, 53, 0, 2.5])

ax2 = ax1.twinx()
ax2.plot(calendar,cases,'go')
ax2.set_ylabel('cases per 10000 people')

plt.xlim([0, 53])
#fig5.suptitle(data.population[region][0], fontsize=20)
fig5.show()
#fig5.savefig('fig5.pdf', bbox_inches='tight')
#exit()
#plt.subplot(2,2,2)
#plt.plot(m,'y--', cases,'g--',ml,'b--',mh,'r--')
#plt.plot(cum_cases,'g--',m_cum,'b--')
################################################################

## Plot Temperature
fig5, ax1 = plt.subplots()

ax1.set_xlabel('time (weeks)')
ax1.set_ylabel('temperature (C)')
ax1.plot(range(0,53),np.array(data.weekly_temperature_array(region,'min')),'b')
ax1.plot(range(0,53),np.array(data.weekly_temperature_array(region,'max')),'r')
ax1.plot(range(0,53),np.array(data.weekly_temperature_array(region,'med')),'y')

ax2 = ax1.twinx()
ax2.plot(calendar,cases,'go')

ax2.set_ylabel('cases per 10000 people')

plt.xlim([0, 53])
#fig5.suptitle(data.population[region][0], fontsize=20)
fig5.show()
input()
#fig5.savefig('fig5temp.pdf', bbox_inches='tight')
