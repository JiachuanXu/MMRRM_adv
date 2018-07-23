	README
===========================
Welcome for any question or suggestion: xjc14@mails.tsinghua.edu.cn

Install:
---------------------------
1. Prerequist: FFTW3, GSL Library, gcc
2. Install:
	Extract the tar file in WORKDIR, then (in terminal)
	$ cd Programs
	$ make test_MMRRM test_powers test_statool test_pdf test_mudecomp test_power_ratio
3. If you need re-compile:
	$ make clean
	$ make test_MMRRM test_powers test_statool test_pdf test_mudecomp test_power_ratio

Quick Start:
---------------------------
1) Edit user changable parameters, flags & definition in 
   	WORKDIR/Programs/header/mmrrm_adv.h
   (You could just turn off all the analysis flags if you only need the output data cube)
   (NOTE: IF YOU DO NOT HAVE Ts DATA, SET BOTH "OPTHIN" AND "HIGHTS"=1, AND SET Ts DIR IN
    WORKDIR/prmt/arguments_mmrrm_adv.dat "NULL")
2) Re-compile the program
3) Edit arguments file 
   	WORKDIR/prmt/args_mmrrm_adv.dat
4) $ cd Programs/
5) $ ./test_MMRRM
6) Now the RSD corrected delta_Tb data cube should be in 
   	WORKDIR/Output_boxes
   (Note: the size of delta_Tb data cube is dim*dim*dim*sizeof(double), so the cube is 
   an double array of shape (dim-x, dim-y, dim-z). )

Source File Configuration:
---------------------------
MMRRM_adv.c: 
	The main program, which mainly operates the MMRRM scheme. There are various of 
	flags in the file "mmrrm_adv.h" in Programs/header, for users to decide whether 
	to turn on some specific function. Some parameters and constants are also defined 
	in "mmrrm_adv.h". The main program would call the following functions (powers, 
	mudecomp, pixelfluxpdf, statool) if needed. That means, you do not need to call 
	these programs deliberately, they are called by test_MMRRM as subroutine.

	Usage: 
		1) Edit the arguments file in ./prmt
		   See "./prmt/args_mmrrm_adv_instruction.txt" for further explanation.
		2) ./test_MMRRM 
	The program would read-in the data cube according to prmt file, then calculate
	RSD corrected data cube and dump them it into ./Output_boxes. After that the main
	program would call subroutines to do some statistic analysis, all the analysis 
	results are dumped into ../Output_PS (by default). For the output format of these 
	results, see following terms.

powers.c: 
	Calculates the auto power spectrum, in k or log(k) interval. To calculate cross ps, 
	some adaptions are necessary. 
	Relies on FFTW3

	Usage(if singly): 
		./test_powers <Tb along x> <Tb along y> <Tb along z>
	<> refers to the file name(WITHOUT PATH, you can edit the path in mmrrm_adv.h) of 
	brightness temperature output, with mesh to mesh applied along x, y or z direction.
	Output format:
	linear: kint k_magnitude(1/Mpc) power(mK^2) deviation
	log:    k_magnitude(1/Mpc) power(mK^2) deviation

mudecomp.c:
	Calculates the fitting ps for \mu^4 term, where \mu refers to the cosine of the 
	angle between 3-D k vector and LoS. (Yi Mao et.al. 2012)
	Relies on GSL(SVD fitting)
	P.S. There's maybe some problems in fitting for log(k) interval, I'll fix it as 
	long as I have time.

	Usage(if singly): ./test_mudecomp <Tb along x> <Tb along y> <Tb along z>
	Output format: 
	linear: k(1/Mpc) mu_4(fitting) dev_mu_4(fitting) chi_square P_mu4(quasi-linear) scatter(quasi-linear)
	log:    k(1/Mpc) mu_4(fitting) dev_mu_4(fitting) chi_square P_mu4(quasi-linear) scatter(quasi-linear)

pixelfluxpdf.c: 
	Calculates the probability distribution function of Tb or other data cube.

	Usage(if singly): ./test_pdf <Tb along x> <Tb along y> <Tb along z> [Output file name]
	[] means optional
	Output format:
	delta_Tb[mK] (counts/total_counts)

statool.c:
	Calculating some statistics: mean, deviation, rms, skewness and kurtosis. in order to
	read data, this file will read an input list, which records the file name of the 
	input data cube. The main program would generate a list in ./Log_files/Output_List 
	by default.

	Usage: ./statool <List name(path also should be included)>

Folders:
---------------------------
Log_files:
	Record the log file for each program
Output_boxes:
	Record the RSD corrected delta_Tb data cube
Output_PS:
	Record the analysis results 
Scripts:
	Record the python scripts used to plot the results as well as the figure
prmt:
	Arguments file for the main program
Programs:
	Source files
