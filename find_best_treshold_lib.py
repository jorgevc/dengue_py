import numpy as np

def is_alert_triggered(m,t,delta,c,stdev):
	if((m[t]-m[t-delta])>c*stdev):
		return True
	else:
		return False
		
def model_alerts(region,data,delta,c,tauVar):
	mosquitoPopulation = mosquitos(region,data)
	weeklyCases = data.cases[region][1:]
	alerts = []
	filtered_alerts = []
	for t in range(len(mosquitoPopulation)):
		stdev = backward_local_std(weeklyCases,t,tauVar)[1]
		if(is_alert_triggered(mosquitoPopulation,t,delta,c,stdev)):
			alerts.append(t)
	for i in range(len(alerts)-1,0,-1):
		if(alerts[i]!=(alerts[i-1]+1)):
			filtered_alerts.append(alerts[i])
	if(len(alerts)>0):
		filtered_alerts.append(alerts[0])	
	return np.array(filtered_alerts[::-1])
	
def outbreaks(cases,delta,e,tauVar):
	outbreaksList = []
	filteredOutbreakList = []
	for t in range(len(cases)):
		if(is_outbreak(cases,t,delta,e,tauVar)):
			outbreaksList.append(t)
	for i in range(len(outbreaksList)-1,0,-1):
		if(outbreaksList[i]!=(outbreaksList[i-1]+1)):
			filteredOutbreakList.append(outbreaksList[i])
	if(len(outbreaksList)>0):		
		filteredOutbreakList.append(outbreaksList[0])
	return np.array(filteredOutbreakList[::-1])
	
def backward_local_std(cases,t,tauVar):
	if(t<tauVar):
		tmpCases = np.append(cases[t-tauVar:],cases[:t])
		stDev = np.std(tmpCases)
		mean = np.mean(tmpCases)
	else:
		stDev = np.std(cases[t-tauVar:t])
		mean = np.mean(cases[t-tauVar:t])
	return np.array([mean,stDev])
	
def mosquitos(region,data):
	import mosquitos_lib as mos
	
	def weekly_mosquitos(m):
		weeklym = np.zeros(53)
		for week in range(53):
			weeklym[week]=m[(week*7.0+3.5)*1.0/dt]
		return weeklym
	
	#dummy initial parameter values
	k = 0.5
	w = 0.05
	epsilon = 0.14*1.05042
	phi = 1.0
	C = 1.0
	pi = 0.12

	#initializing the model for low temperature
	mod = mos.model(k,w,epsilon,phi,pi,C)

	#binding temperature data arrays to models
	print "region %s" % data.temp_min[region][0]
	mod.set_monthly_temperature(data.temp_min[region][1:])
	
	# time units are in days
	dt = 0.01 # time steep (0.01 of a day)
	T = (53)*7*1.0/dt # how many time steeps as function of number of days (weeks*days)
	burnning =  (20)*7*1.0/dt # burnning time
	#burnning evolution
	mBurn = 0.5
	pBurn = 0.3
	for i in range(int(T - burnning) ,int(T)):
		#Set the time steep temperature
		mod.set_model_calendar(i*dt) 
		# Obtain the actual model prediction
		mBurn = mBurn + (mod.k*mod.w*pBurn - mod.epsilon*mBurn)*dt
		pBurn = pBurn + (mod.phi*mBurn*(1.0-pBurn/mod.C)-(mod.pi + mod.w)*pBurn)*dt
		
	#initial conditions
		#defining the numpy needed arrays
	m = np.zeros(T)
	p = np.zeros(T)
	m[0] = mBurn
	p[0] = pBurn
	
	for i in range(1,int(T)):
		#Set the time steep temperature
		mod.set_model_calendar(i*dt) 
		# Obtain the actual model prediction
		m[i] = m[i-1] + (mod.k*mod.w*p[i-1] - mod.epsilon*m[i-1])*dt
		p[i] = p[i-1] + (mod.phi*m[i-1]*(1.0-p[i-1]/mod.C)-(mod.pi + mod.w)*p[i-1])*dt
		
	return weekly_mosquitos(m)

def forward_local_mean(cases,t,tauVar):
	if(t+tauVar > (len(cases)-1)):
		tmpCases = np.append(cases[t:],cases[:(tauVar - (len(cases)-1-t))])
		mean = np.mean(tmpCases)
	else:
		mean = np.mean(cases[t:(t+tauVar)])
	return mean

def is_outbreak(incidence,t,delta,e,tauVar):
	resume = backward_local_std(incidence,t-delta,tauVar)
	pastIncidence = resume[0]
	pastSTD = resume[1]
	futureIncidence = forward_local_mean(incidence,t,tauVar)
	if(futureIncidence > ((1.0 + e.mean)*pastIncidence + e.var*pastSTD)):
		#print "mean future incidence : %f at t=%d" % (futureIncidence, t)
		#print "pastIncidence : %f at t=%d " % (pastIncidence, t)
		#print "pastSTD = %f , e=%f" % (pastSTD , e)
		return True
	else:
		return False
	
