#!/usr/bin/python3.5
import numpy as np
import matplotlib.pyplot as plt
for i in np.arange(0,83,1):
	zrl = np.power(1.02,i)*(6*1.0001+1)-1
	## Begin FLAGS
	num_of_plots = 4 # The number of functions to be plotted
	# file_i = ['file name','color','line style','marker','label']
	file0 = ['../Output_PS/PS_Ratio/Ratio_RSSpace_z%05.2f_dim0256_size0300Mpc_FN-TN' % zrl,'b','-','','Full Treatment over Optical Thin Approx (With Ts)']
	file1 = ['../Output_PS/PS_Ratio/Ratio_RSSpace_z%05.2f_dim0256_size0300Mpc_FH-TH' % zrl,'r','--','','Full Treatment over Optical Thin Approx (Without Ts)']
	files = [file0, file1]
	## End FLAGS
	plt.figure(num = 1,figsize=(8,6)) # Edit figure size and number here
	plt.title('The Ratio of D21^2 for z=%.2f(redshift space)' % zrl,size=20) # Edit figure title here
	plt.xlabel('log(k) [1/Mpc]',size=14) # Edit x-axis label here
	plt.ylabel('log(Delta_21^2) [mK^2]',size=14) # Edit y-axis label here
	#plt.ylim(0.1,1.2)
	for file in files:
		fp = open(file[0],'r')
		k_array = []
		pow_rat = []
		pow_scat = []
		Arg_FLAG = 1
		while Arg_FLAG:
			data = fp.readline().split()
			if not data:
				Arg_FLAG = 0
				break
			k_array.append(data[0])
			pow_rat.append(data[1])
			pow_scat.append(data[2])
		fp.close()
		plt.semilogx(k_array,pow_rat,color=file[1],linestyle=file[2],marker=file[3],label=file[4])
		#plt.errorbar(k_array, pow_spec, pow_spec)
	plt.legend(loc='lower right')
	plt.savefig('images/PS_Ratio/Ratio_RSSpace_log_dim256_size300Mpc_z%05.2f.png' % zrl,format='png') # Edit picture name here
	plt.close()