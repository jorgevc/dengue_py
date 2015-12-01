import numpy as np
import matplotlib.pyplot as plt
import csv
import datos_plot_lib as dpl
import find_best_treshold_lib as fbtl

import matplotlib
matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18)
matplotlib.rc('axes', labelsize=18)

data = dpl.data("DATOS")

region = 12
temp = []
for i in range(0,len(data.cases[region])-1):
	temp.append(data.weekly_temperature(i,region,"min"))

print data.cases[region][0]

#mDiap = fbtl.mosquitos_precipitation_sinTemp(region,data)
mTempPrec= fbtl.mosquitos_precipitation(region,data)
mDiap = fbtl.mosquitos_diapause(region,data)
m = fbtl.mosquitos(region,data)
#m = fbtl.mosquitos_precipitation(region,data)

fig1, ax1 = plt.subplots()
ax1.set_xlabel('time (weeks)')
ax1.set_ylabel('mosquitoes/C')
#ax1.plot(range(0,len(data.cases[region])-1) , data.cases[region][1:], 'b+', range(0,len(data.cases[region])-1), temp, 'r')
g1 = ax1.plot(range(len(data.weekly_precipitation_array(region))), mDiap ,'g', label='Modelo con diapausa')
#g2 = ax1.plot(range(len(data.weekly_precipitation_array(region))), m,'r', label='Modelo con Temperatura s P')
#g3 = ax1.plot(range(len(data.weekly_precipitation_array(region))), mTempPrec,'b', label='Modelo completo')
ax1.legend(loc = 'upper left')
#ax1.legend()
ax2 = ax1.twinx()
ax2.plot(range(0,len(data.cases[region])-1) , data.cases[region][1:], 'go', label='Incidencia')
ax2.set_ylabel('cases per 10000 people')
#ax2.legend()
#ax2.legend(loc='best')
plt.xlim([0, 53])
fig1.suptitle(data.cases[region][0], fontsize=20)

fig1.show()
fig1.savefig('GuerreroDiap.pdf', bbox_inches='tight')

###### fig2
fig2, ax1 = plt.subplots()
ax1.set_xlabel('time (weeks)')
ax1.set_ylabel('mosquitoes/C')
ax1.plot(range(len(data.weekly_precipitation_array(region))), mDiap,'g')
#ax1.plot(range(len(data.weekly_precipitation_array(region))), m,'r')
ax1.legend()
ax2 = ax1.twinx()
ax2.plot(range(len(data.weekly_precipitation_array(region))), data.weekly_precipitation_array(region),'bo', label='Precipitacion')
ax2.set_ylabel('precipitation')
#ax2.legend()
plt.xlim([0, 53])
fig2.suptitle(data.cases[region][0], fontsize=20)
fig2.show()
fig2.savefig('BCSPrecipSTemp2.pdf', bbox_inches='tight')



#### correlaciones
#defining numpy arrays for correlation
cases_per_week = np.zeros(53)
mosquitos_diap = np.zeros(53)
mosquitos_low = np.zeros(53)
mosquitos_comp = np.zeros(53)
shifted_precipitation = np.zeros(53)
corr_diap = np.zeros(53)
corr_precipitation = np.zeros(53)
corr_low = np.zeros(53)
corr_comp = np.zeros(53)

corr_diap_max = 0.0
corr_precipitation_max = 0.0
corr_low_max = 0.0
corr_comp_max = 0.0
#build the point pairs for the correlation

print "Building correlation plots..."

for tau in range(0,53):
	for week in range(1,53+1):
		cases_per_week[week-1] = (1000000.0*data.cases[region][week])/data.population[region][1]
		shifted_precipitation[week-1] = data.weekly_precipitation(week-tau, region)
		mosquitos_diap[week-1] = mDiap[week-tau-1]
		mosquitos_low[week-1] = m[week-tau-1]
		mosquitos_comp[week-1] = mTempPrec[week-tau-1]
	#obtain de correlation for each curve with the corresponding shift 'tau'	
	corr_diap[tau] = np.corrcoef(mosquitos_diap, cases_per_week)[0, 1]
	corr_precipitation[tau] = np.corrcoef(shifted_precipitation,cases_per_week)[0,1]
	corr_low[tau] = np.corrcoef(mosquitos_low, cases_per_week)[0,1]
	corr_comp[tau] = np.corrcoef(mosquitos_comp, cases_per_week)[0,1]
	#find the shift of the max correlation
	if(corr_diap_max < corr_diap[tau]):
		corr_diap_max = corr_diap[tau]
		tau_diap = tau
	if(corr_precipitation_max < corr_precipitation[tau]):
		corr_precipitation_max = corr_precipitation[tau]
		tau_precipitation = tau
	if(corr_low_max < corr_low[tau]):
		corr_low_max = corr_low[tau]
		tau_low = tau
	if(corr_comp_max < corr_comp[tau]):
		corr_comp_max = corr_comp[tau]
		tau_comp = tau
			
print "corr_diap = %f" % corr_diap_max
print "tau_diap = %f" % tau_diap
print "corr_precipitation = %f" % corr_precipitation_max
print "tau_precipitation = %f" % tau_precipitation
print "corr_low = %f" % corr_low_max
print "tau_low = %f" % tau_low
print "corr_com = %f" % corr_comp_max
print "tau_com = %f" % tau_comp

fig3=plt.figure(3)
plt.plot(corr_diap, 'g', label='Modelo con precipitacion s T')
#plt.plot(corr_precipitation,'--', label='Precipitacion')
#plt.plot(corr_comp,'b', label='Modelo completo')
#plt.plot(corr_low,'r', label='Modelo con temperatura s P')

plt.xlim([0, 53])
plt.suptitle(data.cases[region][0], fontsize=20)
plt.ylabel('Correlacion')
plt.xlabel('Retardo (semanas)')
plt.legend()

fig3.savefig('GuerreroCorrsDiap.pdf', bbox_inches='tight')
plt.show()




###### RESULTADOS
#Region Guerrero (12) Diapausa 0.1
#corr_diap = 0.652598
#tau_diap = 12.000000
#corr_precipitation = 0.614457
#tau_precipitation = 7.000000
#corr_low = 0.599848
#tau_low = 5.000000

#Guerrero (12) Precipitacion sin diapausa 
#corr_modprec = 0.618162
#tau_modprec = 2.000000
#corr_precipitation = 0.614457
#tau_precipitation = 7.000000
#corr_low = 0.599848
#tau_low = 5.000000

# Baja California Sur (2) Diapausa 0.05
#corr_diap = 0.804355
#tau_diap = 6.000000
#corr_precipitation = 0.934031
#tau_precipitation = 6.000000
#corr_low = 0.839072
#tau_low = 4.000000

# Baja California Sur (2) Precipitacion (NO diapausa)
#corr_mod_prec = 0.947120
#tau_mod_prec = 4.000000
#corr_precipitation = 0.934031
#tau_precipitation = 6.000000
#corr_low = 0.839072
#tau_low = 4.000000

#Guerrero Comparacion de modelos
#corr_diap = 0.608587
#tau_diap = 2.000000
#corr_precipitation = 0.614457
#tau_precipitation = 7.000000
#corr_low = 0.599848
#tau_low = 5.000000

#Guerrero precipitacion sin temp ni diapausa
#corr_precip_mod = 0.607815
#tau_precip_mod = 2.000000
#corr_precipitation = 0.614457
#tau_precipitation = 7.000000

#Baja California Sur precipitacion sin temp ni diapausa
#corr_precip mod = 0.944761
#tau_precip_mod = 3.000000
#corr_precipitation = 0.934031
#tau_precipitation = 6.000000
