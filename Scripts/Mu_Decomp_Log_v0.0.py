#!/usr/bin/python3.5
import numpy as np
import matplotlib.pyplot as plt
for i in np.arange(0,83,1):
	zrl = np.power(1.02,i)*(6*1.0001+1)-1
	## Begin FLAGS
	# file_i = ['file name','color','line style','marker','label']
	file0 = ['../Output_PS/Mu_decomp/MD_FH_RSSpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'b','-','','Full op-dep & Without $T_s$']
	file1 = ['../Output_PS/Mu_decomp/MD_FN_RSSpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'r','--','','Full op-dep & With $T_s$']
	file2 = ['../Output_PS/2/MD_TH_RSSpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'c','--','','Thin op-dep & Without $T_s$']
	file3 = ['../Output_PS/Mu_decomp/MD_TN_RSSpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'g',':','','Thin op-dep & With $T_s$']
	file4 = ['../Output_PS/2/%06.2f_patchy_mudecomp_from_MMRRM_1x_logk.dat' % zrl,'b','-.','','v0.0']
#	files = [file0, file1, file2, file3]
	files = [file2]
	## End FLAGS
	plt.figure(num = 1,figsize=(8,6)) # Edit figure size and number here
	plt.title('$\mu$ Decomposition: $\Delta_{21}^2$ for z=%.2f(redshift space)' % zrl,size=20) # Edit figure title here
	plt.xlabel('log(k) [1/Mpc]',size=14) # Edit x-axis label here
	plt.ylabel('log($\Delta_{21}^2$) [$mK^2$]',size=14) # Edit y-axis label here
	for file in files:
		fp = open(file[0],'r')
		k_array = []
		pow_spec = []
		pmu4 = []
		Arg_FLAG = 1
		while Arg_FLAG:
			data = fp.readline().split()
			if not data:
				Arg_FLAG = 0
				break
			k_array.append(data[0])
			pow_spec.append(data[1])
			pmu4.append(data[4])
		fp.close()
		plt.semilogx(k_array,pow_spec,color=file[1],linestyle=file[2],marker=file[3],label=file[4])
		plt.semilogx(k_array,pmu4,color='black',linestyle='--',label='quasi-linear adv')
	fp = open(file4[0],'r')
	k_array = []
	pow_spec = []
	pmu4 = []
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		k_array.append(data[1])
		pow_spec.append(data[2])
		pmu4.append(data[3])
	fp.close()
	plt.semilogx(k_array,pow_spec,color=file4[1],linestyle=file4[2],label=file4[4])
	plt.semilogx(k_array,pmu4,color='b',linestyle=':',label='quasi-linear v0.0')
	plt.legend(loc='upper right')
	plt.savefig('images/2/Mu_decomp_RSSpace_log_dim256_size300Mpc_z%05.2f.png' % zrl,format='png') # Edit picture name here
	plt.close()