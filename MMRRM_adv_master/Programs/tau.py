#!/usr/bin/python3.5
import string
import numpy as np
import os
for i in np.arange(0,83,1):
	zrl = np.power(1.02,i)*(6*1.0001+1)-1
	## Begin FLAGS
	# file_i = ['file name','color','line style','marker','label']
	file1 = ['tau_binary/tau_pdf_debug_AproxFN_RESpace_z%05.2f_dim0256_size0300_los1' % zrl]
	file2 = ['tau_binary/tau_pdf_debug_AproxFN_RESpace_z%05.2f_dim0256_size0300_los2' % zrl]
	file3 = ['tau_binary/tau_pdf_debug_AproxFN_RESpace_z%05.2f_dim0256_size0300_los3' % zrl]
	## End FLAGS
	os.system('./test_pdf %s %s %s %s' % (file1[0],file2[0],file3[0],'DBG_tau_pdf_AproxFN_RESpace_z%05.2f_dim0256_size0300.dat' % zrl)) # Influence of (1 - Tcmb / Ts) factor
	