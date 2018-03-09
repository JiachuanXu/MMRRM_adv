#!/usr/bin/python3.5
import numpy as np
import matplotlib.pyplot as plt
for i in np.arange(0,83,1):
	zrl = np.power(1.02,i)*(6*1.0001+1)-1
	## Begin FLAGS
	num_of_plots = 4 # The number of functions to be plotted
	# file_i = ['file name','color','line style','marker','label']
	file0 = ['../Output_PS/Mu_decomp/MD_FH_RSSpace_k_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'b','-','','Full op-dep & Without $T_s$']
	file1 = ['../Output_PS/Mu_decomp/MD_FN_RSSpace_k_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'r','--','','Full op-dep & With $T_s$']
	file2 = ['../Output_PS/dbg_mu_decomp/MD_TH_RSSpace_k_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'c','-.','','Thin op-dep & Without $T_s$']
	file3 = ['../Output_PS/Mu_decomp/MD_TN_RSSpace_k_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'g',':','','Thin op-dep & With $T_s$']
#	files = [file0, file1, file2, file3]
#	files = [file0, file1, file3]
	files = [file2]
	## End FLAGS
	plt.figure(num = 1,figsize=(8,6)) # Edit figure size and number here
	plt.title('$\mu$ Decomposition: $\Delta_{21}^2$ for z=%.2f(redshift space)' % zrl,size=20) # Edit figure title here
	plt.xlabel('k [1/Mpc]',size=14) # Edit x-axis label here
	plt.ylabel('$\Delta_{21}^2$ [$mK^2$]',size=14) # Edit y-axis label here
	for file in files:
		fp = open(file[0],'r')
		k_array = []
		pow_spec = []
		pow_scat = []
		pow_lin  = []
		Arg_FLAG = 1
		while Arg_FLAG:
			data = fp.readline().split()
			if not data:
				Arg_FLAG = 0
				break
			k_array.append(data[0])
			pow_spec.append(data[1])
			pow_scat.append(data[2])
			pow_lin.append(data[4])
		fp.close()
		plt.plot(k_array,pow_spec,color=file[1],linestyle=file[2],marker=file[3],label=file[4])
		plt.plot(k_array,pow_lin,color='black',linestyle=':',label='quasi-linear')
		#plt.errorbar(k_array, pow_spec, pow_scat,color=file[1],linestyle=file[2],marker=file[3],label=file[4])
	plt.legend()
	plt.savefig('images/Mudecomp/Linear/Mu_decomp_TH_RSSpace_lin_dim256_size300Mpc_z%05.2f.png' % zrl,format='png') # Edit picture name here
	plt.close()