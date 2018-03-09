#!/usr/bin/python3.5
import numpy as np
import string
from matplotlib.animation import FuncAnimation
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

## First we have to define the functions we need to plot a single All-in-one figure.
## We divided the figure into three parts: Global Evolution, which indicate where we are; Power Spectrum in log-log coordinate, with the ratio full-treatment over optical thin; the PDF 
## First we'll define the function needed to produce the data points, and wrap in into an interface. Then we creat a figure, and refresh the frame
fig = Figure(figsize = (12, 10))
canvas = FigureCanvas(fig)

def update(frame):
	zrl = np.power(1.02,frame)*(6*1.0001+1)-1
	## Prelude: Document information
	file_GE = ['../Log_files/Global_Evolution/FN_RSSpace_dim0256_size0300_2017_12_6_1:23:36','b','-','Es'] # Filename, color, linestyle, label

	file_PS_RS_TN = ['../Output_PS/SAPS_log/SAPS_TN_RSSpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'#FFA500','--',r'Thin $\tau_\nu$ & With $T_s$']
	file_PS_RS_FN = ['../Output_PS/SAPS_log/SAPS_FN_RSSpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'#FF0000',':',r'Full $\tau_\nu$ & With $T_s$']
	file_PS_RS_TH = ['../Output_PS/SAPS_log/SAPS_TH_RSSpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'#00BFFF','--',r'Thin $\tau_\nu$ & Without $T_s$']
	file_PS_RS_FH = ['../Output_PS/SAPS_log/SAPS_FH_RSSpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'#0000FF',':',r'Full $\tau_\nu$ & Without $T_s$']
	file_PS_RS_seq= [file_PS_RS_TN, file_PS_RS_FN, file_PS_RS_TH, file_PS_RS_FH]
	file_RA_RS_wt = ['../Output_PS/PS_Ratio/Ratio_RSSpace_z%05.2f_dim0256_size0300Mpc_FN-TN' % zrl,'b',':',r'With $T_s$']
	file_RA_RS_wot= ['../Output_PS/PS_Ratio/Ratio_RSSpace_z%05.2f_dim0256_size0300Mpc_FH-TH' % zrl,'r','--',r'Without $T_s$']
	file_PDF_RS_TN= ['../Output_PS/Pixelflux_PDF/PDF_TN_RSSpace_z%05.2f_size0300_dim0256_bin100' % zrl,'#FFA500','--',r'Thin $\tau_\nu$ & With $T_s$']
	file_PDF_RS_FN= ['../Output_PS/Pixelflux_PDF/PDF_FN_RSSpace_z%05.2f_size0300_dim0256_bin100' % zrl,'#FF0000',':',r'Full $\tau_\nu$ & With $T_s$']
	file_PDF_RS_TH= ['../Output_PS/Pixelflux_PDF/PDF_TH_RSSpace_z%05.2f_size0300_dim0256_bin100' % zrl,'#00BFFF','--',r'Thin $\tau_\nu$ & Without $T_s$']
	file_PDF_RS_FH= ['../Output_PS/Pixelflux_PDF/PDF_FH_RSSpace_z%05.2f_size0300_dim0256_bin100' % zrl,'#0000FF',':',r'Full $\tau_\nu$ & Without $T_s$']
	file_PDF_RS_seq=[file_PDF_RS_TN, file_PDF_RS_FN, file_PDF_RS_TH, file_PDF_RS_FH]

	file_PS_RE_TN = ['../Output_PS/SAPS_log/SAPS_TN_RESpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'#FFA500','--',r'Thin $\tau_\nu$ & With $T_s$']
	file_PS_RE_FN = ['../Output_PS/SAPS_log/SAPS_FN_RESpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'#FF0000',':',r'Full $\tau_\nu$ & With $T_s$']
	file_PS_RE_TH = ['../Output_PS/SAPS_log/SAPS_TH_RESpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'#00BFFF','--',r'Thin $\tau_\nu$ & Without $T_s$']
	file_PS_RE_FH = ['../Output_PS/SAPS_log/SAPS_FH_RESpace_logk_z%05.2f_dim0256_size0300Mpc.dat' % zrl,'#0000FF',':',r'Full $\tau_\nu$ & Without $T_s$']
	file_PS_RE_seq= [file_PS_RE_TN, file_PS_RE_FN, file_PS_RE_TH, file_PS_RE_FH]
	file_RA_RE_wt = ['../Output_PS/PS_Ratio/Ratio_RESpace_z%05.2f_dim0256_size0300Mpc_FN-TN' % zrl,'b',':',r'With $T_s$']
	file_RA_RE_wot= ['../Output_PS/PS_Ratio/Ratio_RESpace_z%05.2f_dim0256_size0300Mpc_FH-TH' % zrl,'r','--',r'Without $T_s$']
	file_PDF_RE_TN= ['../Output_PS/Pixelflux_PDF/PDF_TN_RESpace_z%05.2f_size0300_dim0256_bin100' % zrl,'#FFA500','--',r'Thin $\tau_\nu$ & With $T_s$']
	file_PDF_RE_FN= ['../Output_PS/Pixelflux_PDF/PDF_FN_RESpace_z%05.2f_size0300_dim0256_bin100' % zrl,'#FF0000',':',r'Full $\tau_\nu$ & With $T_s$']
	file_PDF_RE_TH= ['../Output_PS/Pixelflux_PDF/PDF_TH_RESpace_z%05.2f_size0300_dim0256_bin100' % zrl,'#00BFFF','--',r'Thin $\tau_\nu$ & Without $T_s$']
	file_PDF_RE_FH= ['../Output_PS/Pixelflux_PDF/PDF_FH_RESpace_z%05.2f_size0300_dim0256_bin100' % zrl,'#0000FF',':',r'Full $\tau_\nu$ & Without $T_s$']
	file_PDF_RE_seq=[file_PDF_RE_TN, file_PDF_RE_FN, file_PDF_RE_TH, file_PDF_RE_FH]

	##Part 1: Construct the canvas

	axes_GE_opthick = fig.add_axes([0.1,0.15,0.2,0.1])
	axes_GE_xHI     = fig.add_axes([0.1,0.25,0.2,0.1])
	axes_GE_TT      = fig.add_axes([0.1,0.35,0.2,0.25])
	axes_GE_Tb      = fig.add_axes([0.1,0.6,0.2,0.25])
	axes_GE_Tb.set_title('Global Evolution(Redshift Space)')
	axes_GE_opthick.set_xlabel('Redshift')
	axes_GE_opthick.set_ylabel('Fraction')
	axes_GE_opthick.set_ylim(0.00000001*0.9, 0.15)
	axes_GE_opthick.set_yticks([0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1])
	axes_GE_opthick.axhline(y=0.01, c='g', linestyle = '-.', label = '$1\%$')
	axes_GE_opthick.axhline(y=0.00001, c='r', linestyle = ':', label = '$0.001\%$')
	axes_GE_xHI.set_xticks([])
	axes_GE_xHI.set_ylabel('$x_{HI}$')
	axes_GE_TT.set_xticks([])
	axes_GE_TT.set_ylabel('T [mK]')
	axes_GE_Tb.set_xticks([])
	axes_GE_Tb.set_ylabel('T [mK]')
	axes_GE_Tb.axhline(y = 0, c = 'g', linestyle = '-.')

	axes_PDF_RS     = fig.add_axes([0.4,0.1,0.25,0.3])
	axes_PDF_RS.set_title('PDF(Redshift Space)')
	axes_PDF_RS.set_xlabel('$T_b$ [mK]')
	axes_PDF_RS.set_ylabel('Probability')
	axes_PDF_RS.set_ylim(0.00000001*0.9, 1*1.1)
	axes_PDF_RE     = fig.add_axes([0.65,0.1,0.25,0.3])
	axes_PDF_RE.set_title('PDF(Real Space)')
	axes_PDF_RE.set_xlabel('$T_b$ [mK]')
	axes_PDF_RE.set_ylim(0.00000001*0.9, 1*1.1)
	axes_PDF_RE.get_yaxis().set_visible(False)

	axes_RA_RS      = fig.add_axes([0.4,0.5,0.25,0.1])
	axes_RA_RS.set_ylim(0.89,1.12)
	axes_RA_RS.set_xlabel('log(k) [1/Mpc]')
	axes_RA_RS.set_ylabel('Ratio F/T')
	axes_RA_RE      = fig.add_axes([0.65,0.5,0.25,0.1])
	axes_RA_RE.set_ylim(0.89,1.12)
	axes_RA_RE.set_xlabel('log(k) [1/Mpc]')
	axes_RA_RE.get_yaxis().set_visible(False)

	axes_PS_RS      = fig.add_axes([0.4,0.6,0.25,0.3])
	axes_PS_RS.set_title('$\Delta_{21}^2$ (Redshift Space)')
	axes_PS_RS.set_ylim(0.02, 1000*1.1)
	axes_PS_RS.get_xaxis().set_visible(False)
	axes_PS_RS.set_ylabel('$\Delta_{21}^2$ $[mK^2]$')
	axes_PS_RE      = fig.add_axes([0.65,0.6,0.25,0.3])
	axes_PS_RE.set_title('$\Delta_{21}^2$ (Real Space)')
	axes_PS_RE.set_ylim(0.02, 1000*1.1)
	axes_PS_RE.get_xaxis().set_visible(False)
	axes_PS_RE.get_yaxis().set_visible(False)

	## Part 2: Functions Definition
	# Setion 1: Global Evolution (Stay fixed)
	zrl_array = np.array([])
	Ts = np.array([])
	Tcmb = np.array([])
	Tb = np.array([])
	xHIv = np.array([])
	xHIm = np.array([])
	opthick = np.array([])

	GE_Arg_FLAG = 1
	fp = open(file_GE[0],'r')
	while GE_Arg_FLAG:
		data = fp.readline().split()
		if not data:
			GE_Arg_FLAG = 0
			break
		zrl_array = np.append(zrl_array, data[0]) ## redshift array
		Ts = np.append(Ts, data[1]) ## Spin temperature
		Tcmb = np.append(Tcmb, data[2]) ## T CMB
		Tb = np.append(Tb, data[3]) ## T 21cm
		xHIv = np.append(xHIv, data[4]) ## xHI, averaged by volumn
		xHIm = np.append(xHIm, data[5]) ## xHI, averaged by mass
		opthick = np.append(opthick, string.atof(data[6])/16777216) ## number of optical thick cells	
	fp.close()
	GE_opthick_line, = axes_GE_opthick.semilogy(zrl_array, opthick, color = 'b', linestyle = '-', label = 'op-thick cells')
	GE_xHI_line1, = axes_GE_xHI.plot(zrl_array, xHIv, color='r',linestyle='-',label='$<x_{HI}>_{cell}$')
	GE_xHI_line2, = axes_GE_xHI.plot(zrl_array, xHIm, color='b',linestyle='--',label='$<x_{HI}>_{mass}$')
	GE_TT_line1, = axes_GE_TT.semilogy(zrl_array,Ts,color='r',linestyle='-',label='$T_s$ [K]')
	GE_TT_line2, = axes_GE_TT.semilogy(zrl_array,Tcmb,color='b',linestyle='--',label='$T_{CMB}$ [K]')
	GE_Tb_line, = axes_GE_Tb.plot(zrl_array, Tb,color='r',linestyle='-',label='$T_b$ of 21cm line[mK]')
	axes_GE_opthick.legend(loc = 'lower right', fontsize = 8)
	axes_GE_xHI.legend(loc = 'lower right', fontsize = 8)
	axes_GE_TT.legend(loc = 'upper right', fontsize = 8)
	axes_GE_Tb.legend(loc = 'lower right', fontsize = 8)

	# Section 2: Global Evolution (Animated)
	#GE_opthick_vline, = 
	axes_GE_opthick.axvline(x = zrl, c = '#DA70D6', linestyle = ':')
	#GE_xHI_vline, = 
	axes_GE_xHI.axvline(x = zrl, c = '#DA70D6', linestyle = ':')
	#GE_TT_vline, = 
	axes_GE_TT.axvline(x = zrl, c = '#DA70D6', linestyle = ':')
	#GE_Tb_line, = 
	axes_GE_Tb.axvline(x = zrl, c = '#DA70D6', linestyle = ':')

	# Section 3: Probability Distribution Function (Animated)

	for file in file_PDF_RS_seq:
		fp = open(file[0],'r')
		PDF_RS_Tbdif = np.array([])
		PDF_RS_Proba = np.array([])
		PDF_Arg_FLAG = 1
		while PDF_Arg_FLAG:
			data = fp.readline().split()
			if not data:
				PDF_Arg_FLAG = 0
				break
			PDF_RS_Tbdif = np.append(PDF_RS_Tbdif, data[0])
			PDF_RS_Proba = np.append(PDF_RS_Proba, data[1])
		fp.close()
		axes_PDF_RS.semilogy(PDF_RS_Tbdif, PDF_RS_Proba, color=file[1], linestyle=file[2],label = file[3])
	axes_PDF_RS.legend(loc = 'lower right', fontsize = 8)

	for file in file_PDF_RE_seq:
		fp = open(file[0],'r')
		PDF_RE_Tbdif = np.array([])
		PDF_RE_Proba = np.array([])
		PDF_Arg_FLAG = 1
		while PDF_Arg_FLAG:
			data = fp.readline().split()
			if not data:
				PDF_Arg_FLAG = 0
				break
			PDF_RE_Tbdif = np.append(PDF_RE_Tbdif, data[0])
			PDF_RE_Proba = np.append(PDF_RE_Proba, data[1])
		fp.close()
		axes_PDF_RE.semilogy(PDF_RE_Tbdif, PDF_RE_Proba, color=file[1], linestyle=file[2],label=file[3])
	axes_PDF_RE.legend(loc = 'lower right', fontsize = 8)

	# Section 4: Power Spectrum (Animated)
	for file in file_PS_RS_seq:
		fp = open(file[0],'r')
		PS_RS_karray = np.array([])
		PS_RS_PSlog  = np.array([])
		PS_Arg_FLAG = 1
		while PS_Arg_FLAG :
			data = fp.readline().split()
			if not data:
				PS_Arg_FLAG = 0
				break
			PS_RS_karray = np.append(PS_RS_karray, data[0])
			PS_RS_PSlog  = np.append(PS_RS_PSlog, data[1])
		fp.close()
		axes_PS_RS.loglog(PS_RS_karray, PS_RS_PSlog, color=file[1],linestyle=file[2],label=file[3])
	axes_PS_RS.legend(loc = 'lower right', fontsize = 8)

	for file in file_PS_RE_seq:
		fp = open(file[0],'r')
		PS_RE_karray = np.array([])
		PS_RE_PSlog  = np.array([])
		PS_Arg_FLAG = 1
		while PS_Arg_FLAG :
			data = fp.readline().split()
			if not data:
				PS_Arg_FLAG = 0
				break
			PS_RE_karray = np.append(PS_RE_karray, data[0])
			PS_RE_PSlog  = np.append(PS_RE_PSlog, data[1])
		fp.close()
		axes_PS_RE.loglog(PS_RE_karray, PS_RE_PSlog, color=file[1],linestyle=file[2],label=file[3])
	axes_PS_RE.legend(loc = 'lower right', fontsize = 8)

	# Section 5: Ratio (Animated)

	fp = open(file_RA_RS_wot[0],'r')
	RA_RS_karray = np.array([])
	RA_RS_ratio = np.array([])
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		RA_RS_karray = np.append(RA_RS_karray, data[0])
		RA_RS_ratio = np.append(RA_RS_ratio, data[1])
	fp.close()
	axes_RA_RS.semilogx(RA_RS_karray, RA_RS_ratio, color=file_RA_RS_wot[1],linestyle=file_RA_RS_wot[2],label=file_RA_RS_wot[3])

	fp = open(file_RA_RS_wt[0],'r')
	RA_RS_karray = np.array([])
	RA_RS_ratio = np.array([])
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		RA_RS_karray = np.append(RA_RS_karray, data[0])
		RA_RS_ratio = np.append(RA_RS_ratio, data[1])
	fp.close()
	axes_RA_RS.semilogx(RA_RS_karray, RA_RS_ratio, color=file_RA_RS_wt[1],linestyle=file_RA_RS_wt[2],label=file_RA_RS_wt[3])
	axes_RA_RS.legend(loc = 'upper right', fontsize = 8)

	fp = open(file_RA_RE_wot[0],'r')
	RA_RE_karray = np.array([])
	RA_RE_ratio = np.array([])
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		RA_RE_karray = np.append(RA_RE_karray, data[0])
		RA_RE_ratio = np.append(RA_RE_ratio, data[1])
	fp.close()
	axes_RA_RE.semilogx(RA_RE_karray, RA_RE_ratio, color=file_RA_RE_wot[1],linestyle=file_RA_RE_wot[2],label=file_RA_RE_wot[3])

	fp = open(file_RA_RE_wt[0],'r')
	RA_RE_karray = np.array([])
	RA_RE_ratio = np.array([])
	Arg_FLAG = 1
	while Arg_FLAG:
		data = fp.readline().split()
		if not data:
			Arg_FLAG = 0
			break
		RA_RE_karray = np.append(RA_RE_karray, data[0])
		RA_RE_ratio = np.append(RA_RE_ratio, data[1])
	fp.close()
	axes_RA_RE.semilogx(RA_RE_karray, RA_RE_ratio, color=file_RA_RE_wt[1],linestyle=file_RA_RE_wt[2],label=file_RA_RE_wt[3])
	axes_RA_RE.legend(loc = 'upper right', fontsize = 8)

	#fig.savefig('images/AiO/z%5.2f_256_300Mpc.png' % zrl,format='png',dpi = 150)
	#fig.clf()
	return axes_PDF_RE,

animation = FuncAnimation(fig,update, interval = 100, blit = 'Ture', frames = 83)
animation.save('sample.gif', fps=30, dpi=150)














