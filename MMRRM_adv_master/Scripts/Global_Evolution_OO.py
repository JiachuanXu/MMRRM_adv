#!/usr/bin/python3.5
import numpy as np
import string
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
## Begin FLAGS
# file_i = ['file name','color','line style','marker','label']
file = ['../Log_files/Global_Evolution/FN_RSSpace_dim0256_size0300_2017_12_6_1:23:36','b','-','','Es']
## End FLAGS
## Reading data
fp = open(file[0],'r')
zrl_array = []
Ts_array = []
Tcmb_array = []
Tb_array = []
xHIv_array = []
xHIm_array = []
opthick_array = []
Arg_FLAG = 1
while Arg_FLAG:
	data = fp.readline().split()
	if not data:
		Arg_FLAG = 0
		break
	zrl_array.append(data[0]) ## redshift array
	Ts_array.append(data[1]) ## Spin temperature
	Tcmb_array.append(data[2]) ## T CMB
	Tb_array.append(data[3]) ## T 21cm
	xHIv_array.append(data[4]) ## xHI, averaged by volumn
	xHIm_array.append(data[5]) ## xHI, averaged by mass
	opthick_array.append(string.atof(data[6])/16777216) ## number of optical thick cells	
fp.close()
## Ploting
fig = Figure(figsize=(8,10))
canvas = FigureCanvas(fig)
opthick = fig.add_axes([0.15,0.1,0.7,0.15])
xHI = fig.add_axes([0.15,0.25,0.7,0.15])
#xHI.sharex = 'all'
gb = fig.add_axes([0.15,0.4,0.7,0.25])
tb = fig.add_axes([0.15,0.65,0.7,0.25])

line1, = opthick.semilogy(zrl_array,opthick_array,color='b',linestyle='-',label='op-thick cells')
opthick.set_xlabel("redshift")
opthick.set_ylabel("Fraction")
opthick.axhline(y=0.01, c='g', linestyle = '-.', label = '$1\%$ fraction')
opthick.axhline(y=0.00001, c='r', linestyle = ':', label = '$0.001\%$ fraction')
opthick.legend(loc = 'upper right')

line2, = xHI.plot(zrl_array,xHIv_array,color='r',linestyle='-',label='$<x_{HI}>_{cell}$')
line3, = xHI.plot(zrl_array,xHIm_array,color='b',linestyle='--',label='$<x_{HI}>_{mass}$')
xHI.set_ylabel("$x_{HI}$")
xHI.set_xticks([])
xHI.legend(loc = 'lower right')

line3, = gb.semilogy(zrl_array,Ts_array,color='r',linestyle='-',label='$T_s$ [K]')
line4, = gb.semilogy(zrl_array,Tcmb_array,color='b',linestyle='--',label='$T_{CMB}$ [K]')
gb.set_ylabel("T [K]")
gb.set_xticks([])
gb.legend(loc ='upper right')

line5, = tb.plot(zrl_array,Tb_array,color='r',linestyle='-',label='$T_b$ of 21cm line[mK]')
tb.set_ylabel("$T_b$ [mK]")
tb.set_title('Global Evolution(redshift space)',size=20)
tb.set_xticks([])
tb.legend(loc = 'lower right')
tb.axhline(y = 0, c = 'g', linestyle = '-.')

fig.savefig('images/Global_Evolution/FN_RSSpace_256_300Mpc.png',format='png',dpi = 150) # Edit picture name here
