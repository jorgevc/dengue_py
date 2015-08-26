import numpy as np
import csv

class data:
	
	def __init__(self, path):
		self.cases = []
		self.temp_med = []
		self.temp_max = []
		self.temp_min = []
		self.precipitation = []
		self.population = []
		
		with open(path + '/dengue-datos.csv', 'rb') as csvfile:
			spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
			spamreader.next() 
			for tmp_row in spamreader:
				row = [tmp_row[0]]
				for i in range(1, len(tmp_row)-1):
					row.append(float(tmp_row[i]))
				self.cases.append(row)
		
		with open(path + '/temp-med.csv', 'rb') as csvfile:
			spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
			spamreader.next()
			for tmp_row in spamreader:
				row = [tmp_row[0]]
				for i in range(1, len(tmp_row)):
					row.append(float(tmp_row[i]))
				self.temp_med.append(row)
			
		with open(path + '/temp-max.csv', 'rb') as csvfile:
			spamreader = csv.reader(csvfile, delimiter=',', quotechar='|') 
			spamreader.next()
			for tmp_row in spamreader:
				row = [tmp_row[0]]
				for i in range(1, len(tmp_row)):
					row.append(float(tmp_row[i]))
				self.temp_max.append(row)
		
		with open(path + '/temp-min.csv', 'rb') as csvfile:
			spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
			spamreader.next()
			for tmp_row in spamreader:
				row = [tmp_row[0]]
				for i in range(1, len(tmp_row)):
					row.append(float(tmp_row[i]))
				self.temp_min.append(row)
			
		with open(path + '/precipitacion.csv', 'rb') as csvfile:
			spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
			spamreader.next()
			for tmp_row in spamreader:
				row = [tmp_row[0]]
				for i in range(1, len(tmp_row)):
					row.append(float(tmp_row[i]))
				self.precipitation.append(row)
				
		with open(path + '/2014Population.csv', 'rb') as csvfile:
			spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
			for tmp_row in spamreader:
				row = [tmp_row[0]]
				for i in range(1, len(tmp_row)):
					row.append(float(tmp_row[i]))
				self.population.append(row)

	def weekly_temperature(self, week, region , opt2):
		UnitsPerMonth = 4.42
		month = int(week/UnitsPerMonth)
		weekM = week - month*4.42
		month = (month-(month/12)*12)
		
		if(opt2 == 'med'):
			MonthlyTemp = self.temp_med[region][1:]
		if(opt2 == 'max'):
			MonthlyTemp = self.temp_max[region][1:]
		if(opt2 == 'min'):
			MonthlyTemp = self.temp_min[region][1:]
		
		nextMonth = (month + 1) - ((month + 1)/12)*12
		prevMonth = (month - 1)
		if(prevMonth == -1):
			prevMonth = 11
		if(weekM <= 2 ):
			Temperature = MonthlyTemp[prevMonth] + (weekM + 2.42)*(MonthlyTemp[month] - MonthlyTemp[prevMonth])/4.42
		else:
			Temperature = MonthlyTemp[month] + (weekM - 2)*(MonthlyTemp[nextMonth] - MonthlyTemp[month])/4.42
		
		return Temperature
		
	def weekly_precipitation(self, week, region):
		UnitsPerMonth = 4.42
		month = int(week/UnitsPerMonth)
		weekM = week - month*4.42
		month = (month-(month/12)*12)
		
		MonthlyRain = self.precipitation[region][1:]
		
		nextMonth = (month + 1) - ((month + 1)/12)*12
		prevMonth = (month - 1)
		if(prevMonth == -1):
			prevMonth = 11
		if(weekM <= 2 ):
			rain = MonthlyRain[prevMonth] + (weekM + 2.42)*(MonthlyRain[month] - MonthlyRain[prevMonth])/4.42
		else:
			rain = MonthlyRain[month] + (weekM - 2)*(MonthlyRain[nextMonth] - MonthlyRain[month])/4.42
		
		return rain

	def weekly_temperature_array(self,region,opt2):
		temperature = []
		for i in range(1,54):
			temperature.append(self.weekly_temperature(i,region,opt2))
		return temperature
		
	def weekly_precipitation_array(self,region):
		rain = []
		for i in range(1,54):
			rain.append(self.weekly_precipitation(i,region))
		return rain
