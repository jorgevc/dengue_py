import numpy as np
import matplotlib.pyplot as plt
import mosquitos_lib as mos
import datos_plot_lib as dpl
import matplotlib

matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18)
matplotlib.rc('axes', labelsize=18)


k = 0.5
w = 0.05
epsilon = 0.14*1.05042
phi = 1.0
C = 1.0
pi = 0.12
mod = mos.model(k,w,epsilon,phi,pi,C)

L=30
eps = np.zeros(L)
omega = np.zeros(L)
phi = np.zeros(L)
pi = np.zeros(L)
R = np.zeros(L)
temp = np.zeros(L)

for i in range(0,L):
	temp[i] = i + 10
	mod.set_model_temperature(temp[i])
	eps[i] = mod.epsilon
	omega[i] = mod.w
	phi[i] = mod.phi
	pi[i] = mod.pi
	R[i] = 0.5*omega[i]*phi[i]/(eps[i]*(pi[i] + omega[i]))
	

plt.figure(1)
plt.plot(temp,eps)
plt.xlabel('temperature (C)',)
plt.ylabel('Feamale mortality rate')
plt.axis([10, 37, 0, 0.15])
#plt.xlim([10, 37])
#plt.suptitle(, fontsize=20)
plt.savefig('fig1Param.pdf', bbox_inches='tight')


plt.figure(2)
plt.plot(temp,omega)
plt.xlabel('temperature (C)',)
plt.ylabel('Transition rate')
plt.axis([10, 37, 0, 0.14])
plt.savefig('fig2Param.pdf', bbox_inches='tight')

plt.figure(3)
plt.plot(temp,phi)
plt.xlabel('temperature (C)',)
plt.ylabel('Oviposition rate')
plt.axis([10, 37, 0, 10])
plt.savefig('fig3Param.pdf', bbox_inches='tight')

plt.figure(4)
plt.plot(temp,pi)
plt.xlabel('temperature (C)',)
plt.ylabel('Mortality rate (aquatic)')
plt.axis([10, 37, 0, 0.15])
plt.savefig('fig4Param.pdf', bbox_inches='tight')

plt.figure(5)
plt.plot(temp,R)
plt.xlabel('temperature (C)',)
plt.ylabel('Offspring number')
plt.axis([10, 37, 0, 120])
plt.savefig('fig5Param.pdf', bbox_inches='tight')


