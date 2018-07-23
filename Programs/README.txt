The general structure of the code is as following:

MMRRM_adv.c: 
	The main program, which mainly operates the MMRRM scheme. There are various of flags in the file "mmrrm_adv.h" in header/, for users to decide 		whether to turn on some specific function. Some parameters and constants are also defined in "mmrrm_adv.h". The main program would call the 		following functions if needed.

	Usage: ./test_MMRRM (You need to edit the arguments file in ../prmt first, there would be another instruction file)
	Then all the results are generated in ../Output_PS (by default). For the output format, refer to the comments in *.c files.

powers.c: 
	Calculates the auto power spectrum, in k or log(k) interval. To calculate cross ps, some adaptions are necessary. Relies on FFTW3

	Usage(if singly): ./test_powers <Tb along x> <Tb along y> <Tb along z>
	<> refers to the file name(without path, you can edit the path in mmrrm_adv.h) of brightness temperature output, with mesh to mesh applied 		along x, y or z direction.

mudecomp.c:
	Calculates the fitting ps for \mu^4 term, where \mu refers to the cosine of the angle between 3-D k vector and LoS. Relies on GSL(SVD fitting)
	P.S. There's maybe some problems in fitting for log(k) interval, I'll fix it as long as I have time.

	Usage(if singly): ./test_mudecomp <Tb along x> <Tb along y> <Tb along z>

pixelfluxpdf.c: 
	Calculates the probability distribution function of Tb or other data cube.

	Usage(if singly): ./test_pdf <Tb along x> <Tb along y> <Tb along z> [Output file name]

statool.c:
	Calculating some statistics: mean, deviation, rms, skewness and kurtosis. This file will read an input list, in which the file name of the 		input data cube are recorded, to read data. The main program would generate a list in ../Log_files/Output_List by default.

	Usage: ./statool <List name(path also should be included)>

For compilation, execute the following command in terminal

	make test_MMRRM test_powers test_mudecomp test_pdf test_statool (I'll change it into a more elegant method, but now I'm not quite conversant 		with makefile)

And 
	make clean
To delete the exsisting binary executable

Welcome for any question or suggestion: xjc14@mails.tsinghua.edu.cn
