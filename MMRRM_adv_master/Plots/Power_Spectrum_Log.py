#!/usr/bin/python3.5
import numpy as np
import string
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
for i in np.arange(0,83,1):
	zrl = np.power(1.02,i)*(6*1.0001+1)-1
	
	## Redshift Space
	# file_i = ['file name','color','line style','marker','label']
	file_RS0 = ['../Output_PS/SAPS_log/SAPS_TH_RSSpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'b','-','','Thin op-dep & Without Ts']
	file_RS1 = ['../Output_PS/SAPS_log/SAPS_FH_RSSpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'r','--','','Full op-dep & Without Ts']
	file_RS2 = ['../Output_PS/SAPS_log/SAPS_TN_RSSpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'g',':','','Thin op-dep & With Ts']
	file_RS3 = ['../Output_PS/SAPS_log/SAPS_FN_RSSpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'c','-.','','Full op-dep & With Ts']

#	fp = open(file_RS0[0],'r')
#	k_array = np.array([])
#	pow_spec_RS_1 = np.array([])
#	pow_scat_RS_1 = np.array([])
#	Arg_FLAG = 1
#	while Arg_FLAG:
#		data = fp.readline().split()
#		if not data:
#			Arg_FLAG = 0
#			break
#		k_array = np.append(k_array, string.atof(data[0]))
#		pow_spec_RS_1 = np.append(pow_spec_RS_1, string.atof(data[1]))
#		pow_scat_RS_1 = np.append(pow_scat_RS_1, string.atof(data[2]))
#	fp.close()

	fp = open(file_RS1[0],'r')
	pow_spec_RS_2 = np.array([])
	pow_scat_RS_2 = np.array([])
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		pow_spec_RS_2 = np.append(pow_spec_RS_2, string.atof(data[1]))
		pow_scat_RS_2 = np.append(pow_scat_RS_2, string.atof(data[2]))
	fp.close()

	fp = open(file_RS2[0],'r')
	pow_spec_RS_3 = np.array([])
	pow_scat_RS_3 = np.array([])
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		pow_spec_RS_3 = np.append(pow_spec_RS_3, string.atof(data[1]))
		pow_scat_RS_3 = np.append(pow_scat_RS_3, string.atof(data[2]))
	fp.close()

	fp = open(file_RS3[0],'r')
	pow_spec_RS_4 = np.array([])
	pow_scat_RS_4 = np.array([])
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		pow_spec_RS_4 = np.append(pow_spec_RS_4, string.atof(data[1]))
		pow_scat_RS_4 = np.append(pow_scat_RS_4, string.atof(data[2]))
	fp.close()

	## Real Space
	# file_i = ['file name','color','line style','marker','label']
	file_RE0 = ['../Output_PS/SAPS_log/SAPS_TH_RESpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'b','-','','Thin op-dep & Without Ts']
	file_RE1 = ['../Output_PS/SAPS_log/SAPS_FH_RESpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'r','--','','Full op-dep & Without Ts']
	file_RE2 = ['../Output_PS/SAPS_log/SAPS_TN_RESpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'g',':','','Thin op-dep & With Ts']
	file_RE3 = ['../Output_PS/SAPS_log/SAPS_FN_RESpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'c','-.','','Full op-dep & With Ts']

#	fp = open(file_RE0[0],'r')
#	pow_spec_RE_1 = np.array([])
#	pow_scat_RE_1 = np.array([])
#	Arg_FLAG = 1
#	while Arg_FLAG:
#		data = fp.readline().split()
#		if not data:
#			Arg_FLAG = 0
#			break
#		pow_spec_RE_1 = np.append(pow_spec_RE_1, string.atof(data[1]))
#		pow_scat_RE_1 = np.append(pow_scat_RE_1, string.atof(data[2]))
#	fp.close()

	fp = open(file_RE1[0],'r')
	pow_spec_RE_2 = np.array([])
	pow_scat_RE_2 = np.array([])
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		pow_spec_RE_2 = np.append(pow_spec_RE_2, string.atof(data[1]))
		pow_scat_RE_2 = np.append(pow_scat_RE_2, string.atof(data[2]))
	fp.close()

	fp = open(file_RE2[0],'r')
	pow_spec_RE_3 = np.array([])
	pow_scat_RE_3 = np.array([])
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		pow_spec_RE_3 = np.append(pow_spec_RE_3, string.atof(data[1]))
		pow_scat_RE_3 = np.append(pow_scat_RE_3, string.atof(data[2]))
	fp.close()

	fp = open(file_RE3[0],'r')
	pow_spec_RE_4 = np.array([])
	pow_scat_RE_4 = np.array([])
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		pow_spec_RE_4 = np.append(pow_spec_RE_4, string.atof(data[1]))
		pow_scat_RE_4 = np.append(pow_scat_RE_4, string.atof(data[2]))
	fp.close()

	## Ratio data
	file_RS_ra1 = ['../Output_PS/PS_Ratio/Ratio_RSSpace_z%05.2f_dim0256_size0300Mpc_F-T' % zrl,'b','-','',r'F/T (1-$e^{-\tau_\nu}$)']
	file_RS_ra2 = ['../Output_PS/PS_Ratio/Ratio_RSSpace_z%05.2f_dim0256_size0300Mpc_H-N' % zrl,'r','--','',r'H/N (1-$T_{CMB}$/$T_s$)']
	file_RE_ra1 = ['../Output_PS/PS_Ratio/Ratio_RESpace_z%05.2f_dim0256_size0300Mpc_FN-TN' % zrl,'b','-','',r'F/T (1-$e^{-\tau_\nu}$']
	file_RE_ra2 = ['../Output_PS/PS_Ratio/Ratio_RESpace_z%05.2f_dim0256_size0300Mpc_FH-TH' % zrl,'r','--','',r'H/N (1-$T_{CMB}$/$T_s$)']

	fp = open(file_RS_ra1[0],'r')
	ratio_wt_RS = np.array([])
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		ratio_wt_RS = np.append(ratio_wt_RS, string.atof(data[1]))
	fp.close()

	fp = open(file_RS_ra2[0],'r')
	ratio_wot_RS = np.array([])
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		ratio_wot_RS = np.append(ratio_wot_RS, string.atof(data[1]))
	fp.close()

	fp = open(file_RE_ra1[0],'r')
	ratio_wt_RE = np.array([])
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		ratio_wt_RE = np.append(ratio_wt_RE, string.atof(data[1]))
	fp.close()

	fp = open(file_RE_ra2[0],'r')
	ratio_wot_RE = np.array([])
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		ratio_wot_RE = np.append(ratio_wot_RE, string.atof(data[1]))
	fp.close()

	mini_ps = min(pow_spec_RS_2.min(),pow_spec_RS_3.min(),pow_spec_RS_4.min(),pow_spec_RE_1.min(),pow_spec_RE_2.min(),pow_spec_RE_3.min(),pow_spec_RE_4.min())
	maxi_ps = max(pow_spec_RS_2.max(),pow_spec_RS_3.max(),pow_spec_RS_4.max(),pow_spec_RE_1.max(),pow_spec_RE_2.max(),pow_spec_RE_3.max(),pow_spec_RE_4.max())
	mini_ra = min(ratio_wt_RS.min(),ratio_wt_RE.min())
	maxi_ra = max(ratio_wt_RS.max(),ratio_wt_RE.max())
	mini_te = min(ratio_wot_RS.min(),ratio_wot_RE.max())
	maxi_te = max(ratio_wot_RS.max(),ratio_wot_RE.max())

	######## Ploting data
	fig = Figure(figsize=(12,8))
	canvas = FigureCanvas(fig)
	sub_RS_ps = fig.add_axes([0.1,0.3,0.4,0.55])
	sub_RE_ps = fig.add_axes([0.5,0.3,0.4,0.55])
	sub_RS_ra = fig.add_axes([0.1,0.1,0.4,0.1])
	sub_RE_ra = fig.add_axes([0.5,0.1,0.4,0.1])
	sub_RS_te = fig.add_axes([0.1,0.2,0.4,0.1])
	sub_RE_te = fig.add_axes([0.5,0.2,0.4,0.1])

	fig.suptitle('$\Delta_{21}^2$ for z=%.2f' % zrl,size=24)

	sub_RS_ps.set_title('Redshift Space',size = 18)
	sub_RS_ps.set_ylim(mini_ps*0.9,maxi_ps*1.1)
	sub_RS_ps.set_ylabel('$\Delta_{21}^2$ $[mK^2]$',size = 14)
	sub_RS_ps.set_xticks([])
#	ps_RS_TH, = sub_RS_ps.loglog(k_array,pow_spec_RS_1,color=file_RS0[1],linestyle=file_RS0[2],marker=file_RS0[3],label=file_RS0[4])
	ps_RS_FH, = sub_RS_ps.loglog(k_array,pow_spec_RS_2,color=file_RS1[1],linestyle=file_RS1[2],marker=file_RS1[3],label=file_RS1[4])
	ps_RS_TN, = sub_RS_ps.loglog(k_array,pow_spec_RS_3,color=file_RS2[1],linestyle=file_RS2[2],marker=file_RS2[3],label=file_RS2[4])
	ps_RS_FN, = sub_RS_ps.loglog(k_array,pow_spec_RS_4,color=file_RS3[1],linestyle=file_RS3[2],marker=file_RS3[3],label=file_RS3[4])
	sub_RS_ps.legend(loc = 'lower right')
	
	sub_RE_ps.set_title('Real Space',size = 18)
	sub_RE_ps.set_ylim(mini_ps*0.9,maxi_ps*1.1)
	sub_RE_ps.set_xticks(np.array([]))
#	ps_RE_TH, = sub_RE_ps.loglog(k_array,pow_spec_RE_1,color=file_RE0[1],linestyle=file_RE0[2],marker=file_RE0[3],label=file_RE0[4])
	ps_RE_FH, = sub_RE_ps.loglog(k_array,pow_spec_RE_2,color=file_RE1[1],linestyle=file_RE1[2],marker=file_RE1[3],label=file_RE1[4])
	ps_RE_TN, = sub_RE_ps.loglog(k_array,pow_spec_RE_3,color=file_RE2[1],linestyle=file_RE2[2],marker=file_RE2[3],label=file_RE2[4])
	ps_RE_FN, = sub_RE_ps.loglog(k_array,pow_spec_RE_4,color=file_RE3[1],linestyle=file_RE3[2],marker=file_RE3[3],label=file_RE3[4])
	sub_RE_ps.legend(loc = 'lower right')
	sub_RE_ps.get_yaxis().set_visible(False)

	sub_RS_ra.set_ylim(mini_ra*0.9,maxi_ra*1.1)
	sub_RS_ra.set_ylabel('Ratio: F/T', size = 14)
	sub_RS_ra.set_xlabel('log(k) [1/Mpc]', size = 14)
#	ra_RS_H, = sub_RS_ra.semilogx(k_array,ratio_wot_RS,color='r',linestyle='--',label=file_RS_ra2[4])
	ra_RS_N, = sub_RS_ra.semilogx(k_array,ratio_wt_RS,color='b',linestyle=':',label=file_RS_ra1[4])
	sub_RS_ra.legend(loc = 'lower right')

	sub_RE_ra.set_ylim(mini_ra*0.9,maxi_ra*1.1)
	sub_RE_ra.set_yticks(np.array([]))
	sub_RE_ra.set_xlabel('log(k) [1/Mpc]',size = 14)
#	ra_RE_H, = sub_RE_ra.semilogx(k_array,ratio_wot_RE,color='r',linestyle='--',label=file_RE_ra2[4])
	ra_RE_N, = sub_RE_ra.semilogx(k_array,ratio_wt_RE,color='b',linestyle=':',label=file_RE_ra1[4])
	sub_RE_ra.legend(loc = 'lower right')

	sub_RS_te.set_ylim(mini_te*0.9,maxi_te*1.1)
	sub_RS_te.set_ylabel('Ratio: H/N', size = 14)
	sub_RS_te.get_xaxis().set_visible(False)
	sub_RS_te.plot(k_array, ratio_wot_RS, color='r',linestyle='--',label=file_RS_ra2[4])
	sub_RS_te.legend(loc = 'lower right')

	sub_RE_te.set_ylim(mini_te*0.9, maxi_te*1.1)
	sub_RE_te.get_yaxis().set_visible(False)
	sub_RE_te.get_xaxis().set_visible(False)
	sub_RE_te.plot(k_array, ratio_wot_RE, color='r',linestyle='--',label=file_RE_ra2[4])
	sub_RE_te.legend(loc = 'lower right')

	fig.savefig('images/PSLog/PSvsRA_RSvsRE_log_dim256_size300Mpc_z%05.2f.png' % zrl,format='png') # Edit picture name here