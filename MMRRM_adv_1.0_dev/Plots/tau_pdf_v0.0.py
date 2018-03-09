#!/usr/bin/python3.5
import string
import numpy as np
import matplotlib.pyplot as plt
for i in np.arange(0,83,1):
	zrl = np.power(1.02,i)*(6*1.0001+1)-1
	## Begin FLAGS
	# file_i = ['file name','color','line style','marker','label']
	file0 = ['../Output_PS/det_v0.0_dp/pdf_bin_z%05.2f_mmrrm_v0.0.dat' % zrl,'b','--','v0.0']
	file1 = ['../Output_PS/det/DBG_tau_pdf_AproxFN_RSSpace_z%05.2f_dim0256_size0300.dat' % zrl,'r',':','adv']
	## End FLAGS
	plt.figure(num = 1,figsize=(8,6)) # Edit figure size and number here
	plt.title(r'PDF of $|1+\delta_{\partial_r v}|$(z=%.2f)' % zrl,size=20) # Edit figure title here
	plt.xlabel('det',size=14) # Edit x-axis label here
#	plt.xlim(0,0.5);
	plt.ylabel('Probability',size=14) # Edit y-axis label here
	plt.grid()

	fp = open(file0[0],'r')
	Tb_diff = []
	Posib_Den = []
	Arg_FLAG = 1
	fp.readline()
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		Tb_diff.append(data[0])
		Posib_Den.append(data[1])
	fp.close()
	plt.semilogy(Tb_diff,Posib_Den,color=file0[1],linestyle=file0[2],label = file0[3])


	fp = open(file1[0],'r')
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
	plt.semilogy(Tb_diff,Posib_Den,color=file1[1],linestyle=file1[2],label = file1[3])
	plt.xlim(0,2.5)
	plt.ylim(0.00000001,1)
	
	plt.legend(loc = 'upper right')
	plt.savefig('images/PDF/Det_v0.0dp_PDF_RSSpace_dim256_size300_z%05.2f.png' % zrl,format='png') # Edit picture name here
	plt.close()