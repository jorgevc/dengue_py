import numpy as np
from math import sqrt

class model:

	def __init__(self, k, w, epsilon, phi, pi, C):
		self.k = k
		self.w = w
		self.epsilon = epsilon
		self.phi = phi
		self.pi = pi
		self.C = C
		self.Time = 0.0
		self.MonthlyTemp = np.zeros(12)
		
	def set_model_temperature(self,Temperature):
		self.epsilon = self.feamale_mortality_rate(Temperature) 
		self.phi = self.feamale_oviposition_rate(Temperature)
		self.pi = self.aquatic_mortality_rate(Temperature)
		self.w = self.aquatic_transition_rate(Temperature)
		
	def return_model_list(self,Temperature):
		model = {'k' : self.k, 'epsilon' : 0, 'phi' : 0, 'pi' : 0, 'w' : 0, 'R' : 0}
		model['epsilon'] = self.feamale_mortality_rate(Temperature)
		model['phi'] = self.feamale_oviposition_rate(Temperature)
		model['pi'] = self.aquatic_mortality_rate(Temperature)
		model['w'] = self.aquatic_transition_rate(Temperature)
		model['R'] = model['k']*model['w']/(model['pi'] + model['w'])*(model['phi']/model['epsilon'])
		return model
		
	def set_model_calendar(self,FisicalTime):
		self.set_model_temperature(self.calendar_temperature(FisicalTime))
		self.Time = FisicalTime

	def calendar_temperature(self,FisicalTime):
		UnitsPerMonth = 30.91
		month = int(FisicalTime/UnitsPerMonth)
		day = FisicalTime - month*30.91
		month = (month-(month/12)*12)
		prevMonth = (month -1)
		nextMonth= (month + 1) - ((month + 1)/12)*12
		if(prevMonth == -1):
			prevMonth = 11
		if(day < 15):
			Temperature = self.MonthlyTemp[prevMonth] + (day + 15.45)*(self.MonthlyTemp[month] - self.MonthlyTemp[prevMonth])/30.91
		else:
			Temperature = self.MonthlyTemp[month] + (day -15.45)*(self.MonthlyTemp[nextMonth] - self.MonthlyTemp[month])/30.91
		return Temperature

	def set_monthly_temperature(self,temperature_array):
		for i in range(0,12):
			self.MonthlyTemp[i] = float(temperature_array[i])
			#self.MonthlyTemp[i] = 19.5
		#self.MonthlyTemp[0] = 19.1
		#self.MonthlyTemp[1] = 20.0
		#self.MonthlyTemp[2] = 21.2
		#self.MonthlyTemp[3] = 22.9
		#self.MonthlyTemp[4] = 22.5
		#self.MonthlyTemp[5] = 29.0  #29.0
		#self.MonthlyTemp[6] = 22.0 #29.8
		#self.MonthlyTemp[7] = 22.0 #30.3
		#self.MonthlyTemp[8] = 22.0 #29.6
		#self.MonthlyTemp[9] = 22.0 #27.0
		#self.MonthlyTemp[10] = 22.5
		#self.MonthlyTemp[11] = 19.3
	
	def feamale_mortality_rate(self,T):
		a=0.8692
		a1=-0.159
		a2=0.01116
		a3=-0.0003408
		a4=0.000003809
		R=a + a1 * T + a2 * T*T + a3*T*T*T + a4*T*T*T*T
		return R
		
	def feamale_oviposition_rate(self,Temp):
		a=-5.4
		a1=1.8
		a2=-0.2124
		a3=0.01015
		a4=-0.0001515
		R=a + a1 * Temp + a2 * Temp*Temp + a3*Temp*Temp*Temp + a4*Temp*Temp*Temp*Temp
		return R

	def aquatic_mortality_rate(self,Temp):
		a=2.130
		a1=-0.3797
		a2=0.02457
		a3=-0.0006778
		a4=0.000006794
		R=a + a1 * Temp + a2 * Temp*Temp + a3*Temp*Temp*Temp + a4*Temp*Temp*Temp*Temp
		return R

	def aquatic_transition_rate(self,Temp):
		a=0.131
		a1=-0.05723
		a2=0.01164
		a3=-0.001341
		a4=0.00008723
		a5=-0.000003017
		a6=0.00000005153
		a7=-0.000000000342
		R=a + a1 * Temp + a2 * Temp*Temp + a3*Temp*Temp*Temp + a4*Temp*Temp*Temp*Temp + a5*Temp*Temp*Temp*Temp*Temp + a6*Temp*Temp*Temp*Temp*Temp*Temp + a7*Temp*Temp*Temp*Temp*Temp*Temp*Temp
		return R
		
	def R(self):
		R=self.k*self.w/(self.pi + self.w)*(self.phi/self.epsilon)
		return R


	def m0_qs(self):
		m =((self.k*self.w)/self.epsilon)*(self.R()-1.0)/self.R()
		return m

	def p0_qs(self):
		p = (self.R()-1.0)/self.R()
		return p
		
	def p1_qs_h(self,Dm0_qs,Dp0_qs):
		A = (self.phi*self.k*self.w)/(self.epsilon*self.C)
		B = self.phi*(self.m0_qs()/self.C - self.k*self.w/self.epsilon - Dm0_qs/(self.C*self.epsilon)) - self.p0_qs()*self.k*self.w/self.epsilon + (self.phi + self.w)
		C = (self.phi + self.p0_qs())*Dm0_qs/self.epsilon + Dp0_qs
		D = B*B - 4.0*A*C
		if(D>=0):
			p1 = (-B + sqrt(D))/(2.0*A)
		else:
			p1 = (-B - sqrt(-D))/(2.0*A)
		return p1
		
	def Dm0_qs(self,FisicalTime,delta):
		T = delta
		modB = self.return_model_list(self.calendar_temperature(FisicalTime - T))
		modF = self.return_model_list(self.calendar_temperature(FisicalTime))
		MF = ((modF['k']*modF['w'])/modF['epsilon'])*(modF['R']-1.0)/modF['R']
		MB = ((modB['k']*modB['w'])/modB['epsilon'])*(modB['R']-1.0)/modB['R']
		DM = (MF - MB)/T
		return DM
		
	def Dp0_qs(self,FisicalTime):
		T= 28.0
		modB = self.return_model_list(self.calendar_temperature(FisicalTime - T))
		modF = self.return_model_list(self.calendar_temperature(FisicalTime))
		MF = (modF['R']-1.0)/modF['R']
		MB = (modB['R']-1.0)/modB['R']
		DP = (MF - MB)/T
		return DP
		
	def p1_qs(self):
		#return self.p1_qs_h(self.Dm0_qs(self.Time),self.Dp0_qs(self.Time))
		#A = -(self.phi/self.epsilon)*self.Dm0_qs(self.Time)*(1.0-self.p0_qs()) - self.Dp0_qs(self.Time)
		#B = (self.phi/(self.C*self.epsilon))*(self.epsilon*self.m0_qs() - self.Dm0_qs(self.Time)) + (self.pi + self.w)
		A = - self.Dp0_qs(self.Time)/((self.phi*self.m0_qs()/self.C) + (self.pi + self.w))
		return A
		
	def m1_qs(self):
		#m1 = self.p1_qs() - self.Dm0_qs(self.Time)/self.epsilon
		m1 = -self.Dm0_qs(self.Time)/self.epsilon
		return m1
		
