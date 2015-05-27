import numpy as np
import matplotlib.pyplot as plt
import mosquitos_lib as mos
#50*50
dt = 0.001
T = 600.0*1.0/dt

m = np.zeros(T)
p = np.zeros(T)
m0_qs = np.zeros(T)
m1_qs = np.zeros(T)

m[0] = 0
p[0] = 0.3

k = 0.5
w = 0.05
epsilon = 0.14*1.05042
phi = 1.0
C = 1.0
pi = 0.12

mod = mos.model(k,w,epsilon,phi,pi,C)

for i in range(1,int(T)):
	mod.set_model_calendar((i-1)*dt)
	m0_qs[i-1] = mod.p0_qs()
	m1_qs[i-1] = mod.p1_qs()
	m[i] = m[i-1] + (mod.k*mod.w*p[i-1] - mod.epsilon*m[i-1])*dt
	p[i] = p[i-1] + (mod.phi*m[i-1]*(1.0 - p[i-1]/mod.C) - (mod.pi + mod.w)*p[i-1])*dt

plt.plot(m1_qs,'r--',m0_qs,'b--',p,'g--',m/40.0,'y--')
#plt.plot(m,'r--',m_qs,'bs')
plt.show()


