#!/usr/bin/python3.5
import string
import numpy as np
import matplotlib.pyplot as plt
for i in np.arange(0,83,1):
	zrl = np.power(1.02,i)*(6*1.0001+1)-1
	## Begin FLAGS
	# file_i = ['file name','color','line style','marker','label']
	file0 = ['../Output_PS/det/DBG_tau_pdf_AproxFN_RSSpace_z%05.2f_dim0256_size0300.dat' % zrl,'b','-']
	## End FLAGS
	plt.figure(num = 1,figsize=(8,6)) # Edit figure size and number here
	plt.title('Possibility Distribution Function for z=%.2f(redshift space)' % zrl,size=20) # Edit figure title here
	plt.xlabel('det',size=14) # Edit x-axis label here
#	plt.xlim(0,0.5);
	plt.ylabel('Probability',size=14) # Edit y-axis label here
	plt.grid()
	fp = open(file0[0],'r')
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
	plt.semilogy(Tb_diff,Posib_Den,color=file0[1],linestyle=file0[2])
	#plt.errorbar(k_array, pow_spec, pow_scat,color=file[1],linestyle=file[2],marker=file[3],label=file[4])
	plt.savefig('images/PDF/Det_PDF_RSSpace_dim256_size300_z%05.2f.png' % zrl,format='png') # Edit picture name here
	plt.close()