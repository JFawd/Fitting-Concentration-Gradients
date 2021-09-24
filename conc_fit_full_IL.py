"""
File: 	conc_fit_full_IL.py
Author:	Jack Fawdon & Kevin Hurlbutt
Date:	2021-05-
----------------------
This program fits an asymmetric concentration gradient,
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
e = 15631
#f = 2075


def func(x, a, b, c, d, f):
	return  (a)*((b/(3.1415)**0.5)*np.exp(-((((-x+e)/1e6)/b)**2))-((((-x+e)/1e6))*special.erfc(((-x+e)/1e6)/b)))-(d)*((c/(3.1415)**0.5)*np.exp(-(((((x)/1e6))/c)**2))-((((x)/1e6))*special.erfc(((x)/1e6)/c)))+(f/1000)

src_dir = input("Please type the full path to the input directory or input file : ").strip() 

report_list = []
max_list = []
min_list = []
list_dely_plate = []
list_dely_strip = []
a_list = []
b_list = []
c_list = []
d_list = []
f_list = []
a_list_err = []
b_list_err = []
c_list_err = []
d_list_err = []
f_list_err = []
area_plate_list = []
area_strip_list = []
area_plate_err_list = []
area_strip_err_list = []

A = ['4h', '20h', '36h']


k = 0
all_txt_files = sorted(os.listdir(src_dir))
for txt in all_txt_files: 
	print(txt)
	if txt.endswith(".txt") and "max" not in txt and "fit" not in txt and "area" not in txt:
		data = np.loadtxt(src_dir+'/'+txt)
		x=data[:,0]
		ydata=data[:,1]


		gmodel = Model(func)
		params = Parameters()
		params.add('a', value=200, vary=True) #min=100000, max=300000)
		params.add('b', value=0.004, vary=True) #min=1e-3, max=1e-1)
		params.add('c', value=0.004, vary=True, min=1e-6, max=1e-1)
		params.add('d', value=200, vary=True) #min=13000, max=17000)
		params.add('f', value=1000, vary=True)# min=0.9, max=1)
		result = gmodel.fit(ydata, params, x=x, weights=np.sqrt(1.0/ydata))
		color = cm.winter(np.linspace(0,1,5))
		plt.rc('axes', prop_cycle=(cycler('color', color)))
		dely = result.eval_uncertainty(sigma=1)

		m = max((result.best_fit))
		mi = min((result.best_fit))
		max_list.append(m)
		min_list.append(mi)



		plt.plot(x, ydata, 'o')
		plt.plot(x, result.best_fit, label = "Fitted Curve")
		plt.xlabel('Distance / um')
		plt.ylabel('Concentration / $moldm^-3$')

		#ci = result.ci_report()
		#print(ci)
		#b_list.append(result.ci_report())
		#print(dely)

		newfilename = "fit_"+str(txt)
		with open(src_dir+'/'+newfilename,'w') as p:
			for j in range(len(x)):
				p.write(str(x[j]))
				p.write('\t')
				p.write(str(result.best_fit[j]))
				p.write('\n')
				j+1

		a_list.append(result.params['a'].value)
		b_list.append(result.params['b'].value)
		c_list.append(result.params['c'].value)
		d_list.append(result.params['d'].value)
		f_list.append(result.params['f'].value)
		a_list_err.append(result.params['a'].stderr)
		b_list_err.append(result.params['b'].stderr)
		c_list_err.append(result.params['c'].stderr)
		d_list_err.append(result.params['d'].stderr)
		f_list_err.append(result.params['f'].stderr)


		y_eval = gmodel.eval(result.params, x=x)

		y_eval_plate = y_eval[24:]
		x_plate = x[24:]

		y_eval_strip = y_eval[:24]
		x_strip = x[:24]
		
		area_plate = integrate.simps(y_eval_plate, x_plate)
		area_strip = integrate.simps(y_eval_strip, x_strip)

		area_plate_list.append(area_plate)
		area_strip_list.append(area_strip)

		dely = result.eval_uncertainty(sigma=1)

		area_plate_err = integrate.simps(dely[24:], x_plate)
		area_strip_err = integrate.simps(dely[:24], x_strip)

		area_plate_err_list.append(area_plate_err)
		area_strip_err_list.append(area_strip_err)

		plt.fill_between(x, result.best_fit-dely, result.best_fit+dely, color="#ABABAB", label='3-$sigma$ uncertainty band')


		m_dely_plate = dely[47]
		m_dely_strip = dely[0]
		list_dely_plate.append(m_dely_plate)
		list_dely_strip.append(m_dely_strip)

		k+1

with open(src_dir+"/max_min.txt", 'w') as f:
	for i in range(len(max_list)):
		f.write(str(max_list[i]))
		f.write('\t')
		f.write(str(min_list[i]))
		f.write('\t')
		f.write(str(list_dely_plate[i]))
		f.write('\t')
		f.write(str(list_dely_strip[i]))
		f.write('\n')
		i+1

with open(src_dir+"/fitting_parameters.txt", 'w') as f:
	for i in range(len(a_list)):
		f.write(str(a_list[i]))
		f.write('\t')
		f.write(str(b_list[i]))
		f.write('\t')
		f.write(str(c_list[i]))
		f.write('\t')
		f.write(str(d_list[i]))
		f.write('\t')
		f.write(str(f_list[i]))
		f.write('\t')
		f.write(str(a_list_err[i]))
		f.write('\t')
		f.write(str(b_list_err[i]))
		f.write('\t')
		f.write(str(c_list_err[i]))
		f.write('\t')
		f.write(str(d_list_err[i]))
		f.write('\t')
		f.write(str(f_list_err[i]))
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


#ax=plt.gca()
#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)

plt.legend(A)

plt.show()


