import numpy as np
import matplotlib.pyplot as plt
import mosquitos_lib as mos
import datos_plot_lib as dpl
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


region = 2
mod = mos.model(k,w,epsilon,phi,pi,C)
modh = mos.model(k,w,epsilon,phi,pi,C)
modl = mos.model(k,w,epsilon,phi,pi,C)
data = dpl.data("DATOS")
mod.set_monthly_temperature(data.temp_med[region][1:])
modh.set_monthly_temperature(data.temp_max[region][1:])
modl.set_monthly_temperature(data.temp_min[region][1:])

print mod.MonthlyTemp
print modh.MonthlyTemp
print modl.MonthlyTemp


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
	cases[i] = data.cases[region][int(week)+1]/40.0
	cum_cases[i] = cum_cases[i-1] + cases[i]*dt
	m_cum[i] = m_cum[i-1] + ml[i]*dt

plt.subplot(2,2,1)
plt.plot(m,'y--', cases,'g--',ml,'b--',mh,'r--')
#plt.plot(m_cum,'r--', cum_cases/200.0,'g--')
#plt.axis([100, T, 0, 3])	
###############################################################
region = 4
mod = mos.model(k,w,epsilon,phi,pi,C)
modh = mos.model(k,w,epsilon,phi,pi,C)
modl = mos.model(k,w,epsilon,phi,pi,C)
data = dpl.data("DATOS")
mod.set_monthly_temperature(data.temp_med[region][1:])
modh.set_monthly_temperature(data.temp_max[region][1:])
modl.set_monthly_temperature(data.temp_min[region][1:])

print mod.MonthlyTemp
print modh.MonthlyTemp
print modl.MonthlyTemp

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
	cases[i] = data.cases[region][int(week)+1]/10.0
	cum_cases[i] = cum_cases[i-1] + cases[i]*dt
	m_cum[i] = m_cum[i-1] + ml[i]*dt

plt.subplot(2,2,2)
plt.plot(m,'y--', cases,'g--',ml,'b--',mh,'r--')

################################################################
###############################################################
region = 5
mod = mos.model(k,w,epsilon,phi,pi,C)
modh = mos.model(k,w,epsilon,phi,pi,C)
modl = mos.model(k,w,epsilon,phi,pi,C)
data = dpl.data("DATOS")
mod.set_monthly_temperature(data.temp_med[region][1:])
modh.set_monthly_temperature(data.temp_max[region][1:])
modl.set_monthly_temperature(data.temp_min[region][1:])

print mod.MonthlyTemp
print modh.MonthlyTemp
print modl.MonthlyTemp


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
	cases[i] = data.cases[region][int(week)+1]/10.0
	cum_cases[i] = cum_cases[i-1] + cases[i]*dt
	m_cum[i] = m_cum[i-1] + ml[i]*dt

plt.subplot(2,2,3)
plt.plot(m,'y--', cases,'g--',ml,'b--',mh,'r--')

################################################################
###############################################################
region = 6
mod = mos.model(k,w,epsilon,phi,pi,C)
modh = mos.model(k,w,epsilon,phi,pi,C)
modl = mos.model(k,w,epsilon,phi,pi,C)
data = dpl.data("DATOS")
mod.set_monthly_temperature(data.temp_med[region][1:])
modh.set_monthly_temperature(data.temp_max[region][1:])
modl.set_monthly_temperature(data.temp_min[region][1:])

print mod.MonthlyTemp
print modh.MonthlyTemp
print modl.MonthlyTemp


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
	cases[i] = data.cases[region][int(week)+1]/4.0
	cum_cases[i] = cum_cases[i-1] + cases[i]*dt
	m_cum[i] = m_cum[i-1] + ml[i]*dt

plt.subplot(2,2,4)
plt.plot(m,'y--', cases,'g--',ml,'b--',mh,'r--')

################################################################

plt.show()


