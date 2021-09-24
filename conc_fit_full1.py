"""
File: 	conc_fit_full.py
Author:	Jack Fawdon
Date:	2021-07-12
----------------------
This program fits a symmetric concentration gradient,
"""


import os
import matplotlib.pyplot as plt
import matplotlib.pyplot as mpl
import matplotlib.cm as cm
import matplotlib as mpl
from collections import OrderedDict
from lmfit import Model, Parameters, minimize, fit_report
import numpy as np
import math
from scipy import integrate, special
from cycler import cycler

cmaps = OrderedDict()


# distance between electrodes
c = 15000

# initial concentration
d = 1000


def func(x, a, b):
	return  (a)*((b/(3.1415)**0.5)*np.exp(-((((-x+c)/1e6)/b)**2))-((((-x+c)/1e6))*special.erfc(((-x+c)/1e6)/b)))-(a)*((b/(3.1415)**0.5)*np.exp(-(((((x)/1e6))/b)**2))-((((x)/1e6))*special.erfc(((x)/1e6)/b)))+(d/1000)


src_dir = input("Please type the full path to the input directory or input file : ").strip()


#bunch of lists to populate
report_list = []
max_list = []
min_list = []
a_list = []
b_list = []
a_list_err = []
b_list_err = []
area_plate_list = []
area_strip_list = []
area_plate_err_list = []
area_strip_err_list = []
A = ['4h', '12h', '20h', '28h', '36h']



#import file
k = 0
all_txt_files = sorted(os.listdir(src_dir))
for txt in all_txt_files: 
	if txt.endswith(".txt") and "fit" not in txt and "max" not in txt:
		data = np.loadtxt(src_dir+'/'+txt)
		x=data[:,0]
		ydata=data[:,1]

		#lmfit function described above
		gmodel = Model(func)
		params = Parameters()
		params.add('a', value=184000, vary=True) #min=100000, max=300000)
		params.add('b', value=0.004, vary=True, min=1e-6, max=1e-1)
		result = gmodel.fit(ydata, params, x=x, weights=np.sqrt(1.0/ydata))
		color = cm.winter(np.linspace(0,1,5))
		plt.rc('axes', prop_cycle=(cycler('color', color)))


		#taking the max and minimum concentration from the fitting and adding to the list above (used for TDF calculation)
		m = max((result.best_fit))
		mi = min((result.best_fit))
		max_list.append(m)
		min_list.append(mi)

		#plotting the data + fit
		plt.plot(x, ydata, 'o')
		plt.plot(x, result.best_fit, label = "Fitted Curve")
		plt.xlabel('Distance / um')
		plt.ylabel('Concentration / $moldm^-3$')

		#dely = result.eval_uncertainty(sigma=1)
		#plt.fill_between(x, result.best_fit-dely, result.best_fit+dely, color="#ABABAB", label='3-$sigma$ uncertainty band')

		#print(result.fit_report())
	
		#writing the fitted data into a txt file
		newfilename = 'fit'+str(txt)
		with open(src_dir+'/'+newfilename,'w') as p:
			for j in range(len(x)):
				p.write(str(x[j]))
				p.write('\t')
				p.write(str(result.best_fit[j]))
				p.write('\n')
				j+1


		#adding the 'a' (interfacial conc gradient) and 'b'	(diffusion length) and associated errors to the list above	
		a_list.append(result.params['a'].value)
		b_list.append(result.params['b'].value)
		a_list_err.append(result.params['a'].stderr)
		b_list_err.append(result.params['b'].stderr)

		#finding the area under the curve for the plating and stripping side
		y_eval = gmodel.eval(result.params, x=x)

		y_eval_plate = y_eval[25:]
		x_plate = x[25:]

		y_eval_strip = y_eval[:25]
		x_strip = x[:25]
		
		area_plate = integrate.simps(y_eval_plate, x_plate)
		area_strip = integrate.simps(y_eval_strip, x_strip)

		area_plate_list.append(area_plate)
		area_strip_list.append(area_strip)

		#the uncertainty band for the fitting and using this integrl as the error 
		dely = result.eval_uncertainty(sigma=1)

		area_plate_err = integrate.simps(dely[25:], x_plate)
		area_strip_err = integrate.simps(dely[:25], x_strip)

		area_plate_err_list.append(area_plate_err)
		area_strip_err_list.append(area_strip_err)

		k+1

#writing to txt files: max conc; the fitting parameters + error; the areas + error
with open(src_dir+"/max_min.txt", 'w') as f:
	for i in range(len(max_list)):
		f.write(str(max_list[i]))
		f.write('\t')
		i+1
	for i in range(len(max_list)):
		f.write(str(min_list[i]))
		f.write('\n')
		i+1


with open(src_dir+"/fitting_parameters.txt", 'w') as f:
	for i in range(len(a_list)):
		f.write(str(a_list[i]))
		f.write('\t')
		f.write(str(b_list[i]))
		f.write('\t')
		f.write(str(a_list_err[i]))
		f.write('\t')
		f.write(str(b_list_err[i]))
		f.write('\n')
		i+1

with open(src_dir+"/areas.txt", 'w') as f:
	for i in range(len(area_plate_list)):
		f.write(str(area_plate_list[i]))
		f.write('\t')
		f.write(str(area_strip_list[i]))
		f.write('\t')
		f.write(str(area_plate_err_list[i]))
		f.write('\t')
		f.write(str(area_strip_err_list[i]))
		f.write('\n')
		i+1


plt.legend(A)

plt.show()
