#!/usr/bin/python3.5
import string
import numpy as np
import matplotlib.pyplot as plt
for i in np.arange(0,83,1):
	zrl = np.power(1.02,i)*(6*1.0001+1)-1
	## Begin FLAGS
	num_of_plots = 4 # The number of functions to be plotted
	# file_i = ['file name','color','line style','marker','label']
	file0 = ['../Output_PS/Pixelflux_PDF/PDF_FH_RSSpace_z%05.2f_size0300_dim0256_bin100' % zrl,'b','-','','Full op-dep & Without Ts']
	file1 = ['../Output_PS/Pixelflux_PDF/PDF_FN_RSSpace_z%05.2f_size0300_dim0256_bin100' % zrl,'r','--','','Full op-dep & With Ts']
	file2 = ['../Output_PS/Pixelflux_PDF/PDF_TH_RSSpace_z%05.2f_size0300_dim0256_bin100' % zrl,'c','-.','','Thin op-dep & Without Ts']
	file3 = ['../Output_PS/Pixelflux_PDF/PDF_TN_RSSpace_z%05.2f_size0300_dim0256_bin100' % zrl,'g',':','','Thin op-dep & With Ts']
#	files = [file0, file1, file2, file3]
	files = [file0, file1, file3]
	## End FLAGS
	plt.figure(num = 1,figsize=(8,6)) # Edit figure size and number here
	plt.title('Possibility Distribution Function for z=%.2f(redshift space)' % zrl,size=20) # Edit figure title here
	plt.xlabel('Tb [mK]',size=14) # Edit x-axis label here
	plt.ylabel('Probability',size=14) # Edit y-axis label here
	plt.grid()
	for file in files:
		fp = open(file[0],'r')
		Tb_diff = []
		Posib_Den = []
		Arg_FLAG = 1
		while Arg_FLAG:
			data = fp.readline().split()
			if not data:
				Arg_FLAG = 0
				break
			Tb_diff.append(data[0])
			Posib_Den.append(data[1])
		fp.close()
		plt.semilogy(Tb_diff,Posib_Den,color=file[1],linestyle=file[2],marker=file[3],label=file[4])
		#plt.errorbar(k_array, pow_spec, pow_scat,color=file[1],linestyle=file[2],marker=file[3],label=file[4])
	plt.legend(loc='lower left')
	plt.savefig('images/PDF/PDF_for_4sets_RSSpace_dim256_size300_z%05.2f.png' % zrl,format='png') # Edit picture name here
	plt.close()