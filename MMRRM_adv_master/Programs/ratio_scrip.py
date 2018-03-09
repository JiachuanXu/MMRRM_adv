#!/usr/bin/python3.5
import string
import numpy as np
import os
for i in np.arange(0,83,1):
	zrl = np.power(1.02,i)*(6*1.0001+1)-1
	## Begin FLAGS
	# file_i = ['file name','color','line style','marker','label']
	file0 = ['../Output_boxes/Tbx_AproxFH_RESpace_z%05.2f_* ../Output_boxes/Tby_AproxFH_RESpace_z%05.2f_* ../Output_boxes/Tbz_AproxFH_RESpace_z%05.2f_*'  % (zrl,zrl,zrl)]
	file1 = ['../Output_boxes/Tbx_AproxFN_RESpace_z%05.2f_* ../Output_boxes/Tby_AproxFN_RESpace_z%05.2f_* ../Output_boxes/Tbz_AproxFN_RESpace_z%05.2f_*'  % (zrl,zrl,zrl)] # Benchmark, full-treatment of optical depth, and include (1 - Tcmb / Ts) factor
#	file2 = ['../Output_boxes/Tbx_AproxTH_RSSpace_z%05.2f_* ../Output_boxes/Tby_AproxTH_RSSpace_z%05.2f_* ../Output_boxes/Tbz_AproxTH_RSSpace_z%05.2f_*'  % (zrl,zrl,zrl)] # Optical thin, without (1 - Tcmb / Ts) factor. This case does not make sense.
	file3 = ['../Output_boxes/Tbx_AproxTN_RESpace_z%05.2f_* ../Output_boxes/Tby_AproxTN_RESpace_z%05.2f_* ../Output_boxes/Tbz_AproxTN_RESpace_z%05.2f_*'  % (zrl,zrl,zrl)]
	## End FLAGS
#	os.system('./test_power_ratio %s %s %s' % (file0[0],file2[0],'FH-TH'))
#	os.system('./test_power_ratio %s %s %s' % (file1[0],file3[0],'FN-TN'))
	os.system('./test_power_ratio %s %s %s' % (file0[0],file1[0],'H-N')) # Influence of (1 - Tcmb / Ts) factor
	os.system('./test_power_ratio %s %s %s' % (file3[0],file1[0],'T-F')) # Influence of optical thin approximation