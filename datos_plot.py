import numpy as np
import matplotlib.pyplot as plt
import csv
import datos_plot_lib as dpl

data = dpl.data("DATOS")

temp = []
for i in range(0,len(data.cases[2])-1):
	temp.append(data.weekly_temperature(i,2,"med"))

	
plt.plot(range(0,len(data.cases[2])-1) , data.cases[2][1:], 'b', range(0,len(data.cases[2])-1), temp, 'r')
plt.show()

