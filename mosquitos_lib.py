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
		UnitsPerMonth = 30.0
		month = int(FisicalTime/UnitsPerMonth)
		day = FisicalTime - month*30.0
		month = (month-(month/12)*12)
		MonthlyTemp = np.zeros(12)
		MonthlyTemp[0]=13.9
		MonthlyTemp[1]=15.1
		MonthlyTemp[2]=17.1
		MonthlyTemp[3]=19.0
		MonthlyTemp[4]=19.8
		MonthlyTemp[5]=19.4
		MonthlyTemp[6]=18.4
		MonthlyTemp[7]=18.4
		MonthlyTemp[8]=18.2
		MonthlyTemp[9]=17.3
		MonthlyTemp[10]=15.8
		MonthlyTemp[11]=14.5
		nextMonth= (month + 1) - ((month + 1)/12)*12
		Temperature = MonthlyTemp[month] + day*(MonthlyTemp[nextMonth] - MonthlyTemp[month])/30.0
		return Temperature

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
		R=a + a1 * Temp + a2 * Temp*Temp + a3*Temp*Temp*Temp + a4*Temp*Temp*Temp*Temp
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
		
	def Dm0_qs(self,FisicalTime):
		T = 28.0
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
		return self.p1_qs_h(self.Dm0_qs(self.Time),self.Dp0_qs(self.Time))
		
	def m1_qs(self):
		m1 = self.p1_qs() - self.Dm0_qs(self.Time)/self.epsilon
		return m1
