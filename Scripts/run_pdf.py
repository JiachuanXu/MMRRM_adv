#!/usr/bin/python3.5
import string
import numpy as np
import os
for i in np.arange(0,83,1):
	zrl = np.power(1.02,i)*(6*1.0001+1)-1
	## Begin FLAGS
	# file_i = ['file name','color','line style','marker','label']
	file0 = ['../Output_boxes/Tbx_FH_z%05.2f_* ../Output_boxes/Tby_FH_z%05.2f_* ../Output_boxes/Tbz_FH_z%05.2f_*'  % (zrl,zrl,zrl),'PDF_FH_z%05.2f_size0300_dim0256_bin100' % zrl]
	file1 = ['../Output_boxes/Tbx_FN_z%05.2f_* ../Output_boxes/Tby_FN_z%05.2f_* ../Output_boxes/Tbz_FN_z%05.2f_*'  % (zrl,zrl,zrl),'PDF_FN_z%05.2f_size0300_dim0256_bin100' % zrl]
	file2 = ['../Output_boxes/Tbx_TH_z%05.2f_* ../Output_boxes/Tby_TH_z%05.2f_* ../Output_boxes/Tbz_TH_z%05.2f_*'  % (zrl,zrl,zrl),'PDF_TH_z%05.2f_size0300_dim0256_bin100' % zrl]
	file3 = ['../Output_boxes/Tbx_TN_z%05.2f_* ../Output_boxes/Tby_TN_z%05.2f_* ../Output_boxes/Tbz_TN_z%05.2f_*'  % (zrl,zrl,zrl),'PDF_TN_z%05.2f_size0300_dim0256_bin100' % zrl]
	## End FLAGS
	os.system('./test_pdf %s %s' % (file0[0],file0[1]))
	os.system('./test_pdf %s %s' % (file1[0],file1[1]))
	os.system('./test_pdf %s %s' % (file2[0],file2[1]))
	os.system('./test_pdf %s %s' % (file3[0],file3[1]))