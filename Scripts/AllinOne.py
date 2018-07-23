# AllinOne.py
# All_in_One.py rebuild
# Compile all the imformation into a single figure
# For discussion only

import numpy as np
import string, os
import matplotlib.pyplot as plt
from matplotlib.figure import Figure 
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from MMRRM_toolkit import Archive, ReadASCIIFolder

#########################
### Define Parameters ###
#########################

DIMENSION = 256 # dimension of datacube
# the 1st order folder
DATA_PDIR = "/Volumes/JX_exFAT/MMRRM_adv/datarepo_MMRRM_adv"
# Global Evolution file
FILE_GEVO = DATA_PDIR + "/Log_files/Global_Evolution/FN_RSSpace_dim0256_\
size0300_2018_7_18_8:44:20"
# Power Spectrum Folder
PS_DIR = DATA_PDIR + "/Output_PS/SAPS_log"
PDF_DIR = DATA_PDIR + "/Output_PS/Pixelflux_PDF"
FIG_DIR = DATA_PDIR + "/Output_PS/Images"

########################
### Define Functions ###
########################

# PlotGlobalEvolution
# plot the global evolution figure
# Input:
#	filename:	global evolution file
#	fig:		target figure
#	axes_loc:	axes location & size, [opthick, xHI, TT, Tb]
#	z_current:	redshift now
#	dim:		dimension of the datacube, default=256
def PlotGlobalEvolution(filename, fig, axes_loc, z_current, dim=256):
	# read file
	# z Ts Tcmb eta delta_Tb xHI_v xHI_m f_opthick
	data = np.genfromtxt(filename, dtype=None, names=True)
	# construct axes
	ax_opthick = fig.add_axes(axes_loc[0])
	ax_xHI = fig.add_axes(axes_loc[1])
	ax_TT = fig.add_axes(axes_loc[2])
	ax_Tb = fig.add_axes(axes_loc[3])
	# plot data
	# optical thick subfigure
	ax_opthick.semilogy(data['z'], data['f_opthick']/float(dim**3.), c='b',
		linestyle='-',label=r'$\tau_\nu$>0.1')
	ax_opthick.set_xlabel('Redshfit')
	ax_opthick.set_ylabel('Fraction')
	ax_opthick.set_ylim(9e-9, 0.15)
	ax_opthick.set_yticks([1e-8,1e-6,1e-4,1e-3,1e-2])
	ax_opthick.axhline(y=0.01, c='g', linestyle = '-.', label = r'$1\%$')
	ax_opthick.axhline(y=0.00001, c='r', linestyle = ':', label = r'$0.001\%$')
	ax_opthick.legend(loc='lower right', fontsize=8)
	# xHI subfigure 
	ax_xHI.plot(data['z'],data['xHI_v'],c='r',linestyle='-.',
		label=r'$\langle x_{HI} \rangle_{volumn}$')
	ax_xHI.plot(data['z'],data['xHI_m'],c='b',linestyle='--',
		label=r'$\langle x_{HI} \rangle_{mass}$')
	ax_xHI.set_xticks([])
	ax_xHI.set_ylabel(r'$x_{HI}$')
	ax_xHI.legend(loc='lower right', fontsize=8)
	# Ts, Tcmb subfigure
	ax_TT.semilogy(data['z'], data['Ts'], c='r', linestyle='-.',label=r'$T_S$')
	ax_TT.semilogy(data['z'], data['Tcmb'], c='b', linestyle='--', 
		label=r'$T_{CMB}$')
	ax_TT.set_xticks([])
	ax_TT.set_ylabel('T [K]')
	ax_TT.legend(loc='upper right', fontsize=8)
	# delta Tb subfigure
	ax_Tb.plot(data['z'], data['delta_Tb'], c='r', linestyle='-')
	ax_Tb.axhline(y=0., c='g', linestyle='-.')
	ax_Tb.set_xticks([])
	ax_Tb.set_ylabel(r'$\delta T_b$ [mK]')

	# add current redshift line
	ax_opthick.axvline(x = z_current, c='#DA70D6', linestyle=':')
	ax_xHI.axvline(x = z_current, c='#DA70D6', linestyle=':')
	ax_TT.axvline(x = z_current, c='#DA70D6', linestyle=':')
	ax_Tb.axvline(x = z_current, c='#DA70D6', linestyle=':')

# PlotPowerSpectrum
# plot power spectrum, spectrum ratio in both real and redshift space
# Inputs:
#	foldername:	name of the folder storing the ps files
#	fig:		target figure
#	axes_loc:	axes location and size, 
#				[ratio_tau_RS,ratio_tau_RE,ratio_ts_RS,ratio_ts_RE,ps_RS,ps_RE]
#	idx:		the index in the z_array
#	zrl:		the target redshift
#	eta:		eta=(1-Ts/Tcmb) at target redshift
def PlotPowerSpectrum(archives, foldername, fig, axes_loc, idx, zrl, eta):
	# read folders
#	archives = ReadASCIIFolder(foldername, "Power Spectrum")
	# check redshift
	for i in [1,2,3,5,6,7]:
		if (zrl - archives[i].z_array[idx])>0.01:
			print("Error: Input argument zrl does not match with archives!")
			exit()
	# read target data
	# power spectrum
	ps_RS_TN = np.genfromtxt(foldername+'/'+archives[1].filenames[idx], names=True, dtype=None)
	ps_RS_FH = np.genfromtxt(foldername+'/'+archives[2].filenames[idx], names=True, dtype=None)
	ps_RS_FN = np.genfromtxt(foldername+'/'+archives[3].filenames[idx], names=True, dtype=None)
	ps_RE_TN = np.genfromtxt(foldername+'/'+archives[5].filenames[idx], names=True, dtype=None)
	ps_RE_FH = np.genfromtxt(foldername+'/'+archives[6].filenames[idx], names=True, dtype=None)
	ps_RE_FN = np.genfromtxt(foldername+'/'+archives[7].filenames[idx], names=True, dtype=None)
	# construct axes
	ax_ratio_tau_RS = fig.add_axes(axes_loc[0])
	ax_ratio_tau_RE = fig.add_axes(axes_loc[1])
	ax_ratio_ts_RS = fig.add_axes(axes_loc[2])
	ax_ratio_ts_RE = fig.add_axes(axes_loc[3])
	ax_PS_RS = fig.add_axes(axes_loc[4])
	ax_PS_RE = fig.add_axes(axes_loc[5])

	# plot data
	
	# plot power spectrum
	# redshift space
	ax_PS_RS.set_title("Redshift Space")
	ax_PS_RS.set_ylabel(r"$\Delta_{21}^2$ [$mK^2$]")
	ax_PS_RS.loglog(ps_RS_TN['k'], ps_RS_TN['ps'],c='b',linestyle=':',
		label=r"$\tau_\nu$ & $(1-T_s/T_{CMB})$")
	ax_PS_RS.loglog(ps_RS_FH['k'],ps_RS_FH['ps']*eta**2.,c='r',linestyle='--',
		label=r"$(1-exp(-\tau_\nu))$ & $\langle (1-T_s/T_{CMB}) \rangle$")
	ax_PS_RS.loglog(ps_RS_FN['k'],ps_RS_FN['ps'],c='black',linestyle='--',
		label=r"$(1-exp(-\tau_\nu))$ & $(1-T_s/T_{CMB})$")
	ax_PS_RS.set_xticks([])
	ax_PS_RS.legend(loc="lower right", fontsize=8)
	# real space
	ax_PS_RE.set_title("Real Space")
	ax_PS_RE.loglog(ps_RE_TN['k'], ps_RE_TN['ps'],c='b',linestyle=':',
		label=r"$\tau_\nu$ & $(1-T_s/T_{CMB})$")
	ax_PS_RE.loglog(ps_RE_FH['k'],ps_RE_FH['ps']*eta**2.,c='r',linestyle='--',
		label=r"$(1-exp(-\tau_\nu))$ & $\langle (1-T_s/T_{CMB}) \rangle$")
	ax_PS_RE.loglog(ps_RE_FN['k'],ps_RE_FN['ps'],c='black',linestyle='--',
		label=r"$(1-exp(-\tau_\nu))$ & $(1-T_s/T_{CMB})$")
	ax_PS_RE.set_xticks([])
	ax_PS_RE.set_ylim(ax_PS_RS.get_ylim())
	ax_PS_RE.get_yaxis().set_visible(False)
	ax_PS_RE.legend(loc="lower right", fontsize=8)
	
	# plot power spectrum ratio: (1-Ts/Tcmb)
	# redshift space
	ax_ratio_ts_RS.semilogx(ps_RS_FN['k'],ps_RS_FN['ps']/(ps_RS_FH['ps']*eta**2.)-1.,
		c='r', linestyle='-.',label=r"Effect of $\eta/\langle \eta \rangle$: $[\Delta_{21,FN}^2/\Delta_{21,FH}^2]-1$")
#	ax_ratio_ts_RS.set_ylabel(r"$\Delta_{21,FN}^2/\Delta_{21,FH}^2$")
	ax_ratio_ts_RS.axhline(0,c='green',linestyle=':')
	ax_ratio_ts_RS.set_xticks([])
	ax_ratio_ts_RS.legend(loc='lower right', fontsize=8)
	# real space
	ax_ratio_ts_RE.semilogx(ps_RE_FN['k'],ps_RE_FN['ps']/(ps_RE_FH['ps']*eta**2.)-1.,
		c='r', linestyle='-.',label=r"Effect of $\eta/\langle \eta \rangle$: $[\Delta_{21,FN}^2/\Delta_{21,FH}^2]-1$")
	ax_ratio_ts_RE.axhline(0,c='green',linestyle=':')
	ax_ratio_ts_RE.set_yticks([])
	ax_ratio_ts_RE.set_xticks([])
	ax_ratio_ts_RE.set_ylim(ax_ratio_ts_RS.get_ylim())
	ax_ratio_ts_RE.legend(loc='lower right', fontsize=8)
	
	# plot power spectrum ratio: (1-exp(-tau))
	# redshift space
	ax_ratio_tau_RS.semilogx(ps_RS_FN['k'],ps_RS_FN['ps']/ps_RS_TN['ps']-1.,c='b',
		linestyle='-.', label=r'Effect of $[1-exp(-\tau_\nu)]/\tau_\nu$: $[\Delta_{21,FN}^2/\Delta_{21,TN}^2]-1$')
#	ax_ratio_tau_RS.set_ylabel(r"$[\Delta_{21,FN}^2/\Delta_{21,TN}^2]-1$")
	ax_ratio_tau_RS.axhline(0.,c='green',linestyle=':')
	ax_ratio_tau_RS.set_xlabel("Log(k) [1/Mpc]")
	ax_ratio_tau_RS.set_ylim([-0.08,0.02])
	ax_ratio_tau_RS.set_yticks([-0.08,-0.04,0.0])
	ax_ratio_tau_RS.set_xticks([0.2,0.4,0.6,0.8,1])
	ax_ratio_tau_RS.legend(loc="lower right", fontsize=8)
	# real space
	ax_ratio_tau_RE.semilogx(ps_RE_FN['k'],ps_RE_FN['ps']/ps_RE_TN['ps']-1.,c='b',
		linestyle='-.', label=r'Effect of $[1-exp(-\tau_\nu)]/\tau_\nu$: $[\Delta_{21,FN}^2/\Delta_{21,TN}^2]-1$')
	ax_ratio_tau_RE.axhline(0.,c='green',linestyle=':')
	ax_ratio_tau_RE.set_yticks([])
	ax_ratio_tau_RE.set_xlabel("Log(k) [1/Mpc]")
	ax_ratio_tau_RE.set_ylim(ax_ratio_tau_RS.get_ylim())
	ax_ratio_tau_RE.set_xticks([0.2,0.4,0.6,0.8,1])
	ax_ratio_tau_RE.legend(loc="lower right", fontsize=8)

# PlotPDF
# Plot the probability distribution function of delta Tb
# Inputs:
#	foldername:		target folder
#	fig:			target figure
#	axes_loc:		axes location & size
#	idx:			index of the target redshift in z_array
#	zrl:			target redshift
#	eta:			eta at target redshift
def PlotPDF(archives, foldername, fig, axes_loc, idx, zrl, eta):
	# Read folders
#	archives = ReadASCIIFolder(foldername, "PDF")
	# check redshift
	for i in [1,2,3,5,6,7]:
		if (zrl - archives[i].z_array[idx])>0.01:
			print("Error: Input argument zrl does not match with archives!")
			exit()
	# read target data
	# redshift space
	pdf_RS_TN = np.genfromtxt(foldername+'/'+archives[1].filenames[idx], names=True, dtype=None)
	pdf_RS_FH = np.genfromtxt(foldername+'/'+archives[2].filenames[idx], names=True, dtype=None)
	pdf_RS_FN = np.genfromtxt(foldername+'/'+archives[3].filenames[idx], names=True, dtype=None)
	# real space
	pdf_RE_TN = np.genfromtxt(foldername+'/'+archives[5].filenames[idx], names=True, dtype=None)
	pdf_RE_FH = np.genfromtxt(foldername+'/'+archives[6].filenames[idx], names=True, dtype=None)
	pdf_RE_FN = np.genfromtxt(foldername+'/'+archives[7].filenames[idx], names=True, dtype=None)
	# construct axes
	ax_pdf_RS = fig.add_axes(axes_loc[0])
	ax_pdf_RE = fig.add_axes(axes_loc[1])

	# plot pdf
	# redshift space
	ax_pdf_RS.set_title("PDF in Redshift Space")
	ax_pdf_RS.semilogy(pdf_RS_TN['delta_Tb'],pdf_RS_TN['counts_fraction'],c='b',
		linestyle=':',label=r"$\tau_\nu$ & $(1-T_s/T_{CMB})$")
	ax_pdf_RS.semilogy(pdf_RS_FH['delta_Tb']*eta,pdf_RS_FH['counts_fraction'],c='r',
		linestyle='--',label=r"$(1-exp(-\tau_\nu))$ & $\langle (1-T_s/T_{CMB}) \rangle$")
	ax_pdf_RS.semilogy(pdf_RS_FN['delta_Tb'],pdf_RS_FN['counts_fraction'],c='black',
		linestyle='--',label=r"$(1-exp(-\tau_\nu))$ & $(1-T_s/T_{CMB})$")
	ax_pdf_RS.set_ylabel("Probability")
	ax_pdf_RS.set_xlabel(r"$\delta T_b$ [mK]")
	ax_pdf_RS.legend(loc="lower right", fontsize=8)
	# real space
	ax_pdf_RE.set_title("PDF in Real Space")
	ax_pdf_RE.semilogy(pdf_RE_TN['delta_Tb'],pdf_RE_TN['counts_fraction'],c='b',
		linestyle=':',label=r"$\tau_\nu$ & $(1-T_s/T_{CMB})$")
	ax_pdf_RE.semilogy(pdf_RE_FH['delta_Tb']*eta,pdf_RE_FH['counts_fraction'],c='r',
		linestyle='--',label=r"$(1-exp(-\tau_\nu))$ & $\langle (1-T_s/T_{CMB}) \rangle$")
	ax_pdf_RE.semilogy(pdf_RE_FN['delta_Tb'],pdf_RE_FN['counts_fraction'],c='black',
		linestyle='--',label=r"$(1-exp(-\tau_\nu))$ & $(1-T_s/T_{CMB})$")
	ax_pdf_RE.set_yticks([])
	ax_pdf_RE.set_xlabel(r"$\delta T_b$ [mK]")
	ax_pdf_RE.set_ylim(ax_pdf_RS.get_ylim())
	ax_pdf_RE.legend(loc="lower right", fontsize=8)

####################
### Main Program ###
####################

### Part I: Construct Canvas

# set the location & size of axes
# Global Evolution
loc_GE_opthick = [0.1,0.15,0.2,0.1] # number of op-thick cells
loc_GE_xHI = [0.1,0.25,0.2,0.1] # evolution of xHI
loc_GE_TT = [0.1,0.35,0.2,0.25] # evolution of Ts Tcmb
loc_GE_Tb = [0.1,0.6,0.2,0.25] # evolution of delta Tb
# Probability Distribution Fuction
loc_PDF_RS = [0.4,0.1,0.25,0.3] # redshift space pdf
loc_PDF_RE = [0.65,0.1,0.25,0.3] # real space pdf
# Power Spectrum (PS)
loc_RATIO_tau_RS = [0.4,0.5,0.25,0.075] # PS ratio in redshift space, tau
loc_RATIO_tau_RE = [0.65,0.5,0.25,0.075] # PS ratio in real space, tau
loc_RATIO_ts_RS = [0.4,0.575,0.25,0.075] # PS ratio in redshift space, (1-Ts/T_cmb)
loc_RATIO_ts_RE = [0.65,0.575,0.25,0.075] # PS ratio in real space, (1-Ts/T_cmb)
loc_PS_RS = [0.4,0.65,0.25,0.25] # PS in redshift space
loc_PS_RE = [0.65,0.65,0.25,0.25] # PS in real space

# Read redshift array & eta array
data = np.genfromtxt(FILE_GEVO, names=True, dtype=None)
redshift_array = data['z']
eta_array = data['eta']
# Read folders
archives_ps = ReadASCIIFolder(PS_DIR, "Power Spectrum")
archives_pdf = ReadASCIIFolder(PDF_DIR, "Power Spectrum")
# loop through redshift
for i,zrl in enumerate(redshift_array):
	fig = Figure(figsize = (12, 10))
	canvas = FigureCanvas(fig)
	fig.suptitle('Snapshot as z=%.2f'%zrl)

	### Part II: Subfigure 1, Global Ebolution

	axes_loc = [loc_GE_opthick, loc_GE_xHI, loc_GE_TT, loc_GE_Tb]
	PlotGlobalEvolution(FILE_GEVO, fig, axes_loc, zrl, DIMENSION)

	### Part III: Subfigure 2, Power spectrum

	axes_loc = [loc_RATIO_tau_RS, loc_RATIO_tau_RE, loc_RATIO_ts_RS, \
	loc_RATIO_ts_RE, loc_PS_RS, loc_PS_RE]
	PlotPowerSpectrum(archives_ps, PS_DIR, fig, axes_loc, i, zrl, eta_array[i])

	### Part IV: Subfigure 3, Probability Distribution Function

	axes_loc = [loc_PDF_RS, loc_PDF_RE]
	PlotPDF(archives_pdf, PDF_DIR, fig, axes_loc, i, zrl, eta_array[i])

	# save figure
	if os.path.exists(FIG_DIR)==False:
		os.system("mkdir %s"%FIG_DIR)
	if os.path.exists(FIG_DIR+"/AiO")==False:
		os.system("mkdir %s"%(FIG_DIR+"/AiO"))
	fig.savefig(FIG_DIR+'/AiO/z%5.2f_256_300Mpc.png'%zrl, format='png', dpi=150)
	fig.clf()
