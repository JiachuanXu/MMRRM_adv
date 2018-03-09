#!/usr/bin/python3.5
import numpy as np
import matplotlib.pyplot as plt
for i in np.arange(0,83,1):
	zrl = np.power(1.02,i)*(6*1.0001+1)-1
## Begin FLAGS
	num_of_plots = 4 # The number of functions to be plotted
# file_i = ['file name','color','line style','marker','label']
	file0 = ['../Output_PS/duibi/SAPS_TH_RSSpace_k_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'b','--','',r'$\tau_\nu$ & $T_s >> T_{CMB}$'] # v0.0
	file1 = ['../Output_PS/SAPS_lin/SAPS_FH_RSSpace_k_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'c','--','',r'$(1-e^{-\tau_\nu})$ & $T_s >> T_{CMB}$']
	file2 = ['../Output_PS/SAPS_lin/SAPS_FN_RSSpace_k_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'r','-.','',r'$(1-e^{-\tau_\nu})$ & ($1-T_{CMB}/T_s$)'] # Standard
	file3 = ['../Output_PS/SAPS_lin/SAPS_TN_RSSpace_k_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'g','-.','',r'$\tau_\nu$ & ($1-T_{CMB}/T_s$)']
	file4 = ['../Output_PS/dbg_v0.0/PSLin/%06.2f_mmrrm256_power.dat' % zrl,'r',':','','v0.0']
#	files = [file0, file1, file2, file3]
#	files = [file0, file1, file2, file3]
	files = [file0]
## End FLAGS
	plt.figure(num = 1,figsize=(8,6)) # Edit figure size and number here
	plt.title(r'$\Delta_{21}^2$ for z=%.2f' % zrl,size=20) # Edit figure title here
	plt.xlabel('k [1/Mpc]',size=14) # Edit x-axis label here
	plt.ylabel(r'$\Delta_{21}^2$ $[mK^2]$',size=14) # Edit y-axis label here
	for file in files:
		fp = open(file[0],'r')
		k_array = []
		pow_spec = []
		pow_scat = []
		Arg_FLAG = 1
		while Arg_FLAG:
			data = fp.readline().split()
			if not data:
				Arg_FLAG = 0
				break
			k_array.append(data[0])
			pow_spec.append(data[1])
			pow_scat.append(data[2])
		fp.close()
		plt.plot(k_array,pow_spec,color=file[1],linestyle=file[2],marker=file[3],label=file[4])
	#plt.errorbar(k_array, pow_spec, pow_scat,color=file[1],linestyle=file[2],marker=file[3],label=file[4])
	fp = open(file4[0],'r')
	k_array = []
	pow_spec = []
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		k_array.append(data[1])
		pow_spec.append(data[2])
	fp.close()
	plt.plot(k_array,pow_spec,color=file4[1],linestyle=file4[2],marker=file4[3],label=file4[4])
#	plt.legend(loc='upper right')
	plt.savefig('images/duibi/PS_benchmark_dim256_size300Mpc_z%05.2f.png' % zrl,format='png') # Edit picture name here
	plt.close()