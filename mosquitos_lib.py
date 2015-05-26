import numpy as np

class model:

	def __init__(self, k, w, epsilon, phi, pi, C):
		self.k = k
		self.w = w
		self.epsilon = epsilon
		self.phi = phi
		self.pi = pi
		self.C = C
		
	def set_model_temperature(self,Temperature):
		self.epsilon = self.feamale_mortality_rate(Temperature) 
		self.phi = self.feamale_oviposition_rate(Temperature)
		self.pi = self.aquatic_mortality_rate(Temperature)
		self.w = self.aquatic_transition_rate(Temperature)
		
	def set_model_calendar(self,FisicalTime):
		self.set_model_temperature(self.calendar_temperature(FisicalTime))

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

