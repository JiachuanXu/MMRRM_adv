# MMRRM_toolkit.py
# This file include the auxilary functions for data analysis
# of the simulation data
# By Jiachuan Xu, 2018.7.11
import numpy as np 
import struct, os, re



# ReadBinary()
# Read the binary datacube
# Inputs:
#	filename:	the dir+filename of target file
#	dtype:		the data type (i=int-4byte, f=float-4byte, d=double-8byte)
#	dim: 		dimension of the datacube, the actual size is dim^3
# Output:
#	data_array:	the tupple object of the datacube
def ReadBinary(filename, dtype='f', dim=256):
	# determine the format string
	total_size = dim*dim*dim
	format_string = "%d%s"%(total_size, dtype)
	# read the file
	with open(filename, mode='rb') as file:
		file_buffer = file.read()
		# check the buffer size
		C_size = len(file_buffer)
		if C_size!=struct.calcsize(format_string):
			print("Error: Reading buffer size is \
uncompatible with binary file size!")
			exit()
		# interpret data as dtype
		data_array = struct.unpack(format_string, file_buffer)
	return data_array

# Archive
# class to keep record of filenames and field name, approximation, e.t.c.
# Initialization:
#	field name:	the name of the datacube, i.e., deltaTb, velocity...
#	space:		observer space, real space or redshift space
#	aprox_1:	approximation 1: optical thin or optical thick
#	aprox_2:	approximation 2: high Ts limit or not
class Archive:
	# initialize
	def __init__(self, _field_name, _space, _aprox_1, _aprox_2):
		self.field_name = _field_name
		self.space = _space
		self.aprox_temprature = _aprox_1
		self.aprox_tau = _aprox_2
		self.terms = 0
		self.filenames = []
		self.z_array = []

	def add_term(self, filename, z):
		self.terms += 1
		self.filenames = np.append(self.filenames, filename)
		self.z_array = np.append(self.z_array, z)

	def sortwithz(self):
		if self.terms>0:
			index = np.argsort(self.z_array)
			self.filenames = self.filenames[index]
			self.z_array = self.z_array[index]

# ReadASCIIFolder
# archive the files in the folder
# Input:
#	folder:		target folder
#	field_name:	name of the datacubes
# Output:
#	archives:	array of Archive class objs
def ReadASCIIFolder(folder, field_name):

	approximations = [0,0,0,0]
	observer_space = [0,0]
	# read flags
	os.system("ls %s/* > temp"%folder)
	with open("temp", mode='r') as file:
		for line in file.readlines():
			filename = line[:-1].split('/')[-1]
			# recognise approximation flag
			aprox = filename.split('_')[1]
			if aprox=='TH':
				approximations[0] = 1
			elif aprox=='TN':
				approximations[1] = 1
			elif aprox=='FH':
				approximations[2] = 1
			elif aprox=='FN':
				approximations[3] = 1
			else:
				print("Error: Can not recognise approximation flag in filename!")
				exit()
			# recognise observer space flag
			space = filename.split('_')[2]
			if space=='RSSpace':
				observer_space[0] = 1
			elif space=='RESpace':
				observer_space[1] = 1
			else:
				print("Error: Can not recognise observer space flag in filename!")
				exit()
	# construct archive list
	num = sum(approximations)*sum(observer_space)
#	print("Approximations x Oberver space = %d x %d"\
#		%(sum(approximations), sum(observer_space)))
	archives = []
	for i in ['RSSpace','RESpace']:
		for j in ['TH','TN','FH','FN']:
			arch = Archive(field_name, i, j[0], j[1])
			archives.append(arch)
	# read filenames into archive list
	with open("temp", mode='r') as file:
		for line in file.readlines():
			filename = line[:-1].split('/')[-1]
			# recognise approximation flag
			aprox = filename.split('_')[1]
			if aprox=='TH':
				id_1 = 0
			elif aprox=='TN':
				id_1 = 1
			elif aprox=='FH':
				id_1 = 2
			elif aprox=='FN':
				id_1 = 3
			else:
				print("Error: Can not recognise approximation flag in filename!")
				exit()
			# recognise observer space flag
			space = filename.split('_')[2]
			if space=='RSSpace':
				id_2 = 0
			elif space=='RESpace':
				id_2 = 1
			else:
				print("Error: Can not recognise observer space flag in filename!")
				exit()
			# add filename
			index = id_1+id_2*4
			z = re.search(r'_z(.*?)_', filename).group(1)
			z = float(z)
			archives[index].add_term(filename, z)
	for archive in archives:
		archive.sortwithz()
	os.system("rm temp")
	return archives
