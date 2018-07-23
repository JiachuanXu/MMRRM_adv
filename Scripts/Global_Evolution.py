#!/usr/bin/python3.5
import numpy as np
import string
import matplotlib.pyplot as plt
## Begin FLAGS
# file_i = ['file name','color','line style','marker','label']
file = ['../Log_files/Global_Evolution/FN_RESpace_dim0256_size0300_2017_12_4_23:32:59','b','-','','Es']
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
	opthick_array.append(string.atof(data[6])+1) ## number of optical thick cells	
fp.close()

fig, ge = plt.subplots(4,1,sharex='all',figsize=(8,16))
fig.subplots_adjust(hspace=0)

#ge[0].axes(0.14,0.65,0.77,0.3)
ge[0].plot(zrl_array,Tb_array,color='r',linestyle='-',label='Tb [mK]')
ge[0].figsize(6.4,4)
#ge[0].ylabel('T [mK]',size=14)
ge[0].legend(loc='lower right')
ge[0].set_title('Global Evolution(real space)',size=20)
ge[0].grid()

#ge[1].axes(0.14,0.35,0.77,0.3)
ge[1].semilogy(zrl_array,Ts_array,color='r',linestyle='-',label='Ts [K]')
ge[1].semilogy(zrl_array,Tcmb_array,color='b',linestyle='--',label='T_CMB [K]')
#ge[1].ylabel('T [K]',size=14)
ge[1].legend(loc='upper right')
ge[1].grid() 

#ge[2].axes(0.14,0.2,0.77,0.15)
ge[2].plot(zrl_array,xHIv_array,color='r',linestyle='-',label='<xHI>_cell')
ge[2].plot(zrl_array,xHIm_array,color='b',linestyle='--',label='<xHI>_mass')
#ge[2].ylabel('xHI',size=14)
ge[2].legend(loc='lower right')

#ge[3].axes(0.14,0.05,0.77,0.15)
ge[3].semilogy(zrl_array,opthick_array,color='b',linestyle='-',label='Number of op-thick cells')
#ge[3].ylabel('Number Counts',size=14)
ge[3].legend(loc='upper right')
ge[3].set_xlabel('redshift',size=14)

plt.savefig('images/Global_Evolution/FN_RESpace_256_300Mpc.png',format='png') # Edit picture name here
plt.close()