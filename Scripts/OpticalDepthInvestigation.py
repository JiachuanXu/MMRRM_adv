# OpticalDepthInvestigation.py
# Investigate the relation between optical depth, Ts and |1+delta_v|
import numpy as np
import matplotlib
matplotlib.use('pdf') # backend
import matplotlib.pyplot as plt
from MMRRM_toolkit import ReadBinary
import os, re, sys, math

#########################
### Define Parameters ###
#########################
DATA_PDIR = "/scratch/xjc14/datarepo_MMRRM_adv/Output_boxes"
DET_DIR = DATA_PDIR + "/det_binary"
TAU_DIR = DATA_PDIR + "/tau_binary"
TS_DIR = "/scratch/xjc14/21cmFAST/Boxes"
DIM = 256
BINS=50
LEVEL=5
IMAGE_DIR = DATA_PDIR + "/Images"

#####################
### Main Programs ###
#####################

### Read folder, we only need FN approximation, RSSpace in this case
# construct z_array, and dictionary
z_array = []
z_dict = []
os.system("ls %s/tau*RSSpace* > temp"%TAU_DIR)
with open("temp", mode='r') as file:
	for line in file.readlines():
		filename = line[:-1].split('/')[-1]
		z = re.search(r'_z(.*?)_', filename).group(1)
		z = float(z)
		if z not in z_array:
			z_array.append(z)
print("Redshift array size: %d"%len(z_array))
for i,zrl in enumerate(z_array):
	z_dict.append((zrl,i))
z_dict = dict(z_dict)

# construct filelist
my_dtype = np.dtype([('z', float), ('tau', '|U200', (3,)), ('det', '|U200', (3,)), 
	('Ts', '|U200')])
filelist = np.zeros(len(z_array), dtype=my_dtype)
# read tau
with open("temp", mode='r') as file:
	for line in file.readlines():
		filename = line.split('/')[-1]
		# read redshift
		z = re.search(r'_z(.*?)_', filename).group(1)
		z = float(z)
		index = z_dict[z]
		# read los
		los = re.search(r"_los(\d)", filename).group(1)
		los = int(los[0])-1
		# record filename
		filelist[index]['tau'][los] = line[:-1]
		filelist[index]['z'] = z
os.system("rm temp")
# read det
os.system("ls %s/det_AproxFN_RSSpace_z* > temp"%DET_DIR)
with open("temp", mode='r') as file:
	for line in file.readlines():
		filename = line.split('/')[-1]
		# read redshift
		z = re.search(r'_z(.*?)_', filename).group(1)
		z = float(z)
		index = z_dict[z]
		# read los
		los = re.search(r'_los(\d)', filename).group(1)
		los = int(los)-1
		# record filename
		filelist[index]['det'][los] = line[:-1]
os.system("rm temp")
# read Ts
os.system("ls %s/Ts_z* > temp"%TS_DIR)
with open("temp", mode='r') as file:
	for line in file.readlines():
		filename = line[:-1].split('/')[-1]
		# read redshift
		z = re.search(r'_z(.*?)_', filename).group(1)
		z = float(z)
		index = z_dict[z]
		# record filename
		filelist[index]['Ts'] = line[:-1]
os.system("rm temp")

### loop through z
for i,zrl in enumerate(z_array):
	print("Processing z = %.2f..."%zrl)
#	data_tau = np.zeros(3*DIM*DIM*DIM, dtype=float)
	data_det = np.zeros(3*DIM*DIM*DIM, dtype=float)
	data_Ts = np.zeros(3*DIM*DIM*DIM, dtype=float)
	data_Ts_temp = ReadBinary(filelist[i]['Ts'], dtype='f',dim=256)
	indexbin_1 = []
	indexbin_2 = []
	indexbin_3 = []
	indexbin_4 = []
	# loop through los
	for los in [0,1,2]:
		# read data, bin according to tau, calculate fraction
		data_tau_temp = ReadBinary(filelist[i]['tau'][los], dtype='d', dim=256)
		data_det_temp = ReadBinary(filelist[i]['det'][los], dtype='d', dim=256)
		for ct in range(DIM*DIM*DIM):
			tau = data_tau_temp[ct]
#			data_tau[los*DIM*DIM*DIM+ct] = tau
			data_det[los*DIM*DIM*DIM+ct] = data_det_temp[ct]
			data_Ts[los*DIM*DIM*DIM+ct] = data_Ts_temp[ct]
			if tau<0.1:
				indexbin_1.append(los*DIM*DIM*DIM+ct)
			elif tau>0.1 and tau<0.5:
				indexbin_2.append(los*DIM*DIM*DIM+ct)
			elif tau>0.5 and tau<1.0:
				indexbin_3.append(los*DIM*DIM*DIM+ct)
			elif tau>1.0:
				indexbin_4.append(los*DIM*DIM*DIM+ct)
	# determine range
	Ts_lim = [min(data_Ts), max(data_Ts)]
	det_lim = [min(data_det), max(data_det)]
	fig, axes_ = plt.subplots(2,2,sharex=True, sharey=True)
	fig.suptitle("z=%.2f"%zrl)
	axes = [axes_[0,0],axes_[0,1],axes_[1,0],axes_[1,1]]
	indexbin = [indexbin_1, indexbin_2, indexbin_3, indexbin_4]
	for j in range(4):
		# bin into histogram
		hist_temp, xedges, yedges = np.histogram2d(data_Ts[indexbin[j]], 
			data_det[indexbin[j]], range=[Ts_lim, det_lim], 
			bins=[np.logspace(math.log(Ts_lim[0],10), math.log(Ts_lim[1],10), 
			BINS+1, endpoint=True),
			np.linspace(det_lim[0], det_lim[1], BINS+1, endpoint=True)])
		percentage = len(indexbin[j])/float(3*DIM*DIM*DIM)*100.
		# plot the data in 4 bins
		# get 2D percentage density
		hist = hist_temp/(3*DIM*DIM*DIM*(math.log(Ts_lim[1],10)-math.log(Ts_lim[0],10))*
			(det_lim[1]-det_lim[0])/BINS/BINS)
		xcoord = np.sqrt(xedges[1:]*xedges[:-1])
		ycoord = (yedges[1:]+yedges[:-1])/2.
		# plot contour
		axes[j].contourf(ycoord, xcoord, hist, LEVEL, cmap='Reds',origin='lower',
			extent=[det_lim[0],det_lim[1],Ts_lim[0],Ts_lim[1]])
		c = axes[j].contour(ycoord,xcoord,hist, LEVEL, origin='lower',colors='black',extent=[det_lim[0],det_lim[1],Ts_lim[0],Ts_lim[1]],linewidths=0.5,linestyles='dashdot')
		plt.clabel(c, inline=True, fontsize=8)
		axes[j].set_yscale('log')
		axes[j].text(0.9,0.1,r'Per: %.2f'%percentage,horizontalalignment='right',
			verticalalignment='center', transform=axes[j].transAxes)
	axes[0].set_ylabel(r"$T_s$ [K]")
	axes[2].set_ylabel(r"$T_s$ [K]")
	axes[2].set_xlabel(r"$|1+\delta_{\partial_r v}|$")
	axes[3].set_xlabel(r"$|1+\delta_{\partial_r v}|$")

	# save figure
	if os.path.exists(IMAGE_DIR)==False:
		os.system("mkdir %s"%IMAGE_DIR)
	if os.path.exists(IMAGE_DIR+"/TauInvestigate")==False:
		os.system("mkdir %s/TauInvestigate"%IMAGE_DIR)
	fig.savefig(IMAGE_DIR+"/TauInvestigate/z%05.2f_256_300Mpc.png"%zrl, 
		format='png',dpi=150)
