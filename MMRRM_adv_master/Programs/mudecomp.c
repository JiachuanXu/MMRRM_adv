#include "header/read_write_box.h"
#include "header/read_args.h"
#include "header/cal_fftw.h"
#include "header/cal_Tb.h"
#include <gsl/gsl_multifit.h>

#define LOG_MAX_VOLUME (unsigned long long)100000
#define mymax(x,y) (x>y? x:y)
/*
	Usage: ./mudecomp <Tb x> <Tb y> <Tb z> <box PATH> <delta> <xHI> <Ts> <dim_rt>
	Output format: 
	linear: k(1/Mpc) mu_4(fitting) dev_mu_4(fitting) chi_square P_mu4(quasi-linear) scatter(quasi-linear)
	log:    k(1/Mpc) mu_4(fitting) dev_mu_4(fitting) chi_square P_mu4(quasi-linear) scatter(quasi-linear)
*/
	FILE *LOG;
	int dim_nbody, dim_rt, mesh2mesh, hights, opthin;
	float box_size, zrl;
void svd_fit_lin(double *power_bin,double *dev_bin, double *mu_bin,int k_start,int k_end,double *fit_coef_bin,double *fit_coef_err,double *chi2,int term);
void svd_fit_log(double *power_bin,double *dev_bin, double *mu_bin,int k_start,int k_end,double *fit_coef_bin,double *fit_coef_err,double *chi2,int term,int nintvl,double dlnk,unsigned long long binsize);
void getname_MD_Lin(char *filename);
void getname_MD_Log(char *filename);

int main(int argc,char *argv[]){
	// Variable declaration
	FILE *MUDECOMP_LIN, *MUDECOMP_LOG;
	char LOG_NAME[M_BOXNAME], mudecomp_lin[M_BOXNAME], mudecomp_log[M_BOXNAME], cmd[M_CMD];
	double *Tb, *D2_decomp, *sqr_decomp, *mu_decomp;
	double *c2_lin, *dev_lin, *chi2_lin;
	double *c2_log, *dev_log, *chi2_log, *karray;
	double k_step, dlnk, krange, nintvl;
	unsigned long long *count_decomp;
	unsigned long long binsize;
	int idir, mid;
	fftw_complex *power;
	time_t start_time;
	struct tm *local;
	start_time=time(NULL);
	local=localtime(&start_time);

	if(argc!=9){
		printf("USAGE: ./test_mudecomp <Tb x> <Tb y> <Tb z> <box PATH> <delta> <xHI> <Ts> <dim_rt>\n");
		exit(EXIT_FAILURE);
	}
	system("mkdir ../Log_files");
	system("mkdir ../Log_files/Mu_decomp");
	sprintf(LOG_NAME,"../Log_files/Mu_decomp/mudecomp_log_file_%d_%d_%d_%d:%d:%d",local->tm_year+YEAR_START,local->tm_mon+MON_START,local->tm_mday,local->tm_hour,local->tm_min,local->tm_sec);
	LOG=fopen(LOG_NAME,"w");
	read_title(argv[1],&zrl,&box_size,&dim_nbody,&opthin,&hights,&mesh2mesh);
	dim_rt = atoi(argv[8]);
	printf("****** Calling Mu Decompositon: z = %.2f\n",zrl);
	mid = dim_nbody/2;//dim is even
	krange = mid - ORDER/2;// For linear case
	binsize = (unsigned long long)(mid+1+ORDER/2)*(mid-ORDER/2)/2;
	k_step = 2.0*PI/(double)box_size;
	nintvl = (int)(log((double)(mid-1)/kmincut)/log((kmincut+1.0)/kmincut));
	dlnk = (double)(log((mid-1)/kmincut)/((double)nintvl-1.0));
	if(MUDECINK){
		c2_lin = dmyalloc(krange);	// coefficient of the 4 order polynomical, linear
		dev_lin = dmyalloc(krange);	// deviation of c2_lin
		chi2_lin = dmyalloc(krange);	// chi square
	}
	if(MUDECINLOGK){
		c2_log = dmyalloc(nintvl);	// coefficient of the 4 order polynomical, log
		dev_log = dmyalloc(nintvl);	// deviation of c2_log
		chi2_log = dmyalloc(nintvl);	// chi square
		karray = dmyalloc(nintvl);	// k [1/Mpc]
	}
	// Azimuthal averaged data input for fitting
	D2_decomp = dmyalloc(binsize);
	sqr_decomp = dmyalloc(binsize);
	mu_decomp = dmyalloc(binsize);
	count_decomp = ullmyalloc(binsize);
	/******************************** 	SVD fitting: produce azimuthal averaged data (D_21, mu, deviation) 	********************************/
	for(idir=1;idir<4;idir++){
		printf("Looping through LoS %d\n",idir);
		Tb = dmyalloc(box_vol(dim_nbody));
		power = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*kbox_vol(dim_nbody));
		dbox_read(Tb,MMRRM_BOX_OP,argv[idir],dim_nbody);
		calculate_power(Tb,Tb,dim_nbody,power);
		azimuthal_sorting(power, count_decomp, D2_decomp, sqr_decomp, dim_nbody, mid, idir);
		fftw_free(power); power = NULL;
		free(Tb); Tb = NULL;
	}
	azimuthal_average(count_decomp, D2_decomp, sqr_decomp, mu_decomp, mid);
	free(count_decomp); count_decomp = NULL;
	// Now the standard azimuthal averaged data have been produced, refer data point of (kint, mu) = (l,m) with bin_index(l,m)
	if(MUDECINK)
		svd_fit_lin(D2_decomp,sqr_decomp,mu_decomp,ORDER/2,mid,c2_lin,dev_lin,chi2_lin,4);
	if(MUDECINLOGK)
		svd_fit_log(D2_decomp,sqr_decomp,mu_decomp,mymax(ORDER/2, kmincut),mid,c2_log,dev_log,chi2_log,4,nintvl,dlnk,binsize);
	free(D2_decomp); D2_decomp = NULL;
	free(sqr_decomp); sqr_decomp = NULL;
	free(mu_decomp); mu_decomp = NULL;
	/******************************** Quasi-linear decomposition: P_mu4 = ((delta Tb hat)(eta_mean))^2*D_HH ********************************/
	float xHIm_ave;
	double eta_bar, Tb_h;
	double *power_spec_Pmu4_lin, *power_spec_Pmu4_log, *err_Pmu4_lin, *err_Pmu4_log, *k_array_lin, *k_array_log;
	unsigned long long *num_Pmu4_lin, *num_Pmu4_log;

	float *delta, *xHI, *Ts;
	double *d_delta;
	fftw_complex *power_Pmu4;
	delta = fmyalloc(box_vol(dim_nbody));
	d_delta = dmyalloc(box_vol(dim_nbody));
	xHI = fmyalloc(box_vol(dim_rt));
	Ts = fmyalloc(box_vol(dim_nbody));
	fbox_read(delta, argv[4],argv[5],dim_nbody);
	fbox_read(xHI,argv[4],argv[6],dim_rt);
	xHIm_ave = Mass_Ave(delta, xHI, dim_nbody, dim_rt);
	Tb_h = Tb_hat(zrl, xHIm_ave);
	for(unsigned long long ct = 0; ct < box_vol(dim_nbody); ct++)
		d_delta[ct] = (double)delta[ct];
	free(delta); delta = NULL;
	free(xHI); xHI = NULL;
	fbox_read(Ts,argv[4],argv[7],dim_nbody);
	eta_bar = Eta_Ave(Ts, zrl, dim_nbody);
	free(Ts); Ts = NULL;

	power_Pmu4 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*kbox_vol(dim_nbody));
	calculate_power(d_delta,d_delta,dim_nbody,power_Pmu4);
	free(d_delta); d_delta = NULL;
	if(MUDECINK){
		power_spec_Pmu4_lin = dmyalloc(mid - 1);
		err_Pmu4_lin = dmyalloc(mid - 1);
		k_array_lin = dmyalloc(mid - 1);
		num_Pmu4_lin = ullmyalloc(mid - 1);
		spherical_sorting_linear(power_Pmu4, num_Pmu4_lin, power_spec_Pmu4_lin, err_Pmu4_lin, dim_nbody);
		spherical_average_linear(num_Pmu4_lin, power_spec_Pmu4_lin, err_Pmu4_lin, k_array_lin, dim_nbody, mid, k_step);
		free(num_Pmu4_lin); num_Pmu4_lin = NULL;
	}
	if(MUDECINLOGK){
		power_spec_Pmu4_log = dmyalloc(nintvl);
		err_Pmu4_log = dmyalloc(nintvl);
		k_array_log = dmyalloc(nintvl);
		num_Pmu4_log = ullmyalloc(nintvl);
		spherical_sorting_log(power_Pmu4, num_Pmu4_log, power_spec_Pmu4_log, err_Pmu4_log, dim_nbody, nintvl, dlnk);
		spherical_average_log(num_Pmu4_log, power_spec_Pmu4_log, err_Pmu4_log, k_array_log, dim_nbody, nintvl, dlnk, k_step);
		free(num_Pmu4_log); num_Pmu4_log = NULL;
	}
	fftw_free(power_Pmu4); power_Pmu4 = NULL;
	
	/******************************** 						  Print the result  						 ********************************/
	sprintf(cmd,"mkdir %s",MMRRM_PS_OP);
	system(cmd);
	sprintf(cmd,"mkdir %s/Mu_decomp",MMRRM_PS_OP);
	system(cmd);
	if(MUDECINK){
		getname_MD_Lin(mudecomp_lin);
		MUDECOMP_LIN = fopen(mudecomp_lin,"w");
		fprintf(MUDECOMP_LIN,"# k\tmu4_fit\tmu4_fit_err\tchi2\tPmu4\tscatter\n");
		for(int ct=ORDER/2;ct<mid;ct++){
			c2_lin[ct - ORDER/2] *= 4.0*PI*pow((double)ct,3.0)/pow((double)dim_nbody,6.0);
			dev_lin[ct - ORDER/2] *= 4.0*PI*pow((double)ct,3.0)/pow((double)dim_nbody,6.0);
			// k(1/Mpc) mu_4(fitting) dev_mu_4(fitting) chi_square mu_4(quasi-linear) dev_mu_4(quasi-linear)
			fprintf(MUDECOMP_LIN,"%E\t%E\t%E\t%E\t%E\t%E\n",k_array_lin[ct - 1],c2_lin[ct - ORDER/2],dev_lin[ct - ORDER/2],chi2_lin[ct - ORDER/2], (Tb_h * eta_bar)*(Tb_h * eta_bar)*power_spec_Pmu4_lin[ct - 1], (Tb_h * eta_bar)*(Tb_h * eta_bar)*err_Pmu4_lin[ct - 1]);
		}
		fclose(MUDECOMP_LIN);
		free(c2_lin); c2_lin = NULL;
		free(dev_lin); dev_lin = NULL;
		free(chi2_lin); chi2_lin = NULL;
		free(power_spec_Pmu4_lin); power_spec_Pmu4_lin = NULL;
		free(err_Pmu4_lin); err_Pmu4_lin = NULL;
		free(k_array_lin); k_array_lin = NULL;
	}
	if(MUDECINLOGK){
		getname_MD_Log(mudecomp_log);
		MUDECOMP_LOG = fopen(mudecomp_log,"w");
		fprintf(MUDECOMP_LOG,"# k\tmu4_fit\tmu4_fit_err\tchi2\tPmu4\tscatter\n");
		for(int ct_log=0;ct_log<nintvl;ct_log++){
			karray[ct_log] = exp(dlnk*(double)ct_log+log((double)mymax(kmincut,ORDER/2)));
			c2_log[ct_log] *= 4.0*PI*pow((double)karray[ct_log],3.0)/pow((double)dim_nbody,6.0);
			dev_log[ct_log] *= 4.0*PI*pow((double)karray[ct_log],3.0)/pow((double)dim_nbody,6.0);
			// k(1/Mpc) mu_4(fitting) dev_mu_4(fitting) chi_square mu_4(quasi-linear) dev_mu_4(quasi-linear)
			fprintf(MUDECOMP_LOG,"%E\t%E\t%E\t%E\t%E\t%E\n",karray[ct_log]*k_step,c2_log[ct_log],dev_log[ct_log],chi2_log[ct_log], (Tb_h * eta_bar)*(Tb_h * eta_bar)*power_spec_Pmu4_log[ct_log], (Tb_h * eta_bar)*(Tb_h * eta_bar)*err_Pmu4_log[ct_log]);
		}
		fclose(MUDECOMP_LOG);
		free(c2_log); c2_log = NULL;
		free(dev_log); dev_log = NULL;
		free(chi2_log); chi2_log = NULL;
		free(power_spec_Pmu4_log); power_spec_Pmu4_log = NULL;
		free(err_Pmu4_log); err_Pmu4_log = NULL;
		free(karray); karray = NULL;
		free(k_array_log); k_array_log = NULL;
	}
	printf("Normal ending...\n");
	fprintf(LOG,"Normal ending...\n");
	fclose(LOG);
	return 0;
}

void svd_fit_lin(double *power_bin,double *dev_bin, double *mu_bin, int k_start,int k_end,double *fit_coef_bin,double *fit_coef_err,double *chi2,int term){
	int l, nfit = (ORDER/2 + 1);
	#pragma omp parallel num_threads(nthreads) default(none) shared(k_start,k_end,nfit,mu_bin,power_bin,dev_bin,term,fit_coef_bin,fit_coef_err,chi2) private(l)
	{
		#pragma omp for
		for(l = k_start;l < k_end;l++){
			gsl_matrix *X,*cov;
			gsl_vector *y,*w,*c;
			int ndata = l+1;
			double chisq;
			X = gsl_matrix_alloc(ndata,nfit);
			y = gsl_vector_alloc(ndata);
			w = gsl_vector_alloc(ndata);
			c = gsl_vector_alloc(nfit);
			cov = gsl_matrix_alloc(nfit,nfit);
			// Setting vectors and matrix
			for(int m = 0;m <= l;m++){
				for(int pm = 0;pm <= ORDER/2;pm++){
					gsl_matrix_set(X,m,pm,pow(mu_bin[bin_index(l,m)],2.0*pm));
				}
				gsl_vector_set(y,m,power_bin[bin_index(l,m)]);
				gsl_vector_set(w,m,1.0/(dev_bin[bin_index(l,m)]*dev_bin[bin_index(l,m)]));
			}
			gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(ndata,nfit);
			gsl_multifit_wlinear(X,w,y,c,cov,&chisq,work);
			//#pragma omp atomic
			fit_coef_bin[l - k_start] = gsl_vector_get(c,term/2);
			//#pragma omp atomic
			fit_coef_err[l - k_start] = sqrt(gsl_matrix_get(cov,2,2));
			//#pragma omp atomic
			chi2[l - k_start] = chisq;
			gsl_multifit_linear_free(work); work = NULL;
			gsl_matrix_free(X); X = NULL;
			gsl_vector_free(y); y = NULL;
			gsl_vector_free(w); w = NULL;
			gsl_vector_free(c); c = NULL;
			gsl_matrix_free(cov); cov = NULL;
		}
	}
}

void svd_fit_log(double *power_bin,double *dev_bin, double *mu_bin,int k_start,int k_end,double *fit_coef_bin,double *fit_coef_err,double *chi2,int term,int nintvl,double dlnk,unsigned long long binsize){
	int nfit = ORDER/2+1,l,l_log;
	unsigned long long ct_stake[nintvl + 1];
	int ndata[nintvl];
	for(int i=0;i<nintvl;i++)
		ndata[i] = 0;
	// Setting ndata array
	for(l = k_start;l < k_end;l++){
		int kint_log = (int)(log((double)l/(double)kmincut)/dlnk + 0.5);
		if(kint_log >= 0||kint_log < nintvl)
			ndata[kint_log] += (l + 1);
		else
			printf("*** mudecomp.c: ERROR: mudecomp.c: kint_log out of range!\n");
	}
	ct_stake[0] = (ORDER/2 + 1 + k_start)*(k_start - ORDER/2)/2;
	for(int i=1;i < nintvl + 1; i++){
			ct_stake[i] = ct_stake[i-1] + ndata[i-1];
	}
	#pragma omp parallel num_threads(nthreads) default(none) shared(nintvl, nfit,term,fit_coef_bin,fit_coef_err,chi2, mu_bin, power_bin, dev_bin, ndata, ct_stake) private(l_log)
	{
		#pragma omp for
		for(l_log=0;l_log < nintvl;l_log++){
			gsl_matrix *X,*cov;
			gsl_vector *y,*w,*c;
			double chisq;
			X = gsl_matrix_alloc(ndata[l_log],nfit);
			y = gsl_vector_alloc(ndata[l_log]);
			w = gsl_vector_alloc(ndata[l_log]);
			c = gsl_vector_alloc(nfit);
			cov = gsl_matrix_alloc(nfit,nfit);

			for(unsigned long long ct = ct_stake[l_log];ct < ct_stake[l_log + 1];ct++){
				unsigned long long ct_dif = ct - ct_stake[l_log];
				for(int pm = 0;pm <= ORDER/2;pm++){
					if(ct_dif < ndata[l_log])
						gsl_matrix_set(X,ct_dif,pm,pow(mu_bin[ct],2.0*pm));
					else
						printf("*** mudecomp.c: Error: mudecomp.c: Index out of range\nct:%llu\n",ct);
				}
				gsl_vector_set(y,ct_dif,power_bin[ct]);
				gsl_vector_set(w,ct_dif,1.0/(dev_bin[ct] * dev_bin[ct]));
			}
			gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(ndata[l_log],nfit);
			gsl_multifit_wlinear(X,w,y,c,cov,&chisq,work);
			fit_coef_bin[l_log] += gsl_vector_get(c,term/2);
			fit_coef_err[l_log] += sqrt(gsl_matrix_get(cov,2,2));
			chi2[l_log] += chisq;
			gsl_multifit_linear_free(work); work = NULL;
			gsl_matrix_free(X); X = NULL;
			gsl_vector_free(y); y = NULL;
			gsl_vector_free(w); w = NULL;
			gsl_vector_free(c); c = NULL;
			gsl_matrix_free(cov); cov = NULL;
		}
	}
}

void getname_MD_Lin(char *filename){
	char label_Aprox_O, label_Aprox_T;
	char label_rr[8] = {'\0'};
	(hights)? (label_Aprox_T = 'H'):(label_Aprox_T = 'N');
	(opthin)? (label_Aprox_O = 'T'):(label_Aprox_O = 'F');
	(mesh2mesh)? strcpy(label_rr, "RSSpace"):strcpy(label_rr,"RESpace");
	sprintf(filename,"%s/Mu_decomp/MD_%c%c_%s_k_z%05.2f_dim%04d_size%04.0fMpc.dat",MMRRM_PS_OP,label_Aprox_O,label_Aprox_T,label_rr,zrl,dim_nbody,box_size);
}
void getname_MD_Log(char *filename){
	char label_Aprox_O, label_Aprox_T;
	char label_rr[8] = {'\0'};
	(hights)? (label_Aprox_T = 'H'):(label_Aprox_T = 'N');
	(opthin)? (label_Aprox_O = 'T'):(label_Aprox_O = 'F');
	(mesh2mesh)? strcpy(label_rr, "RSSpace"):strcpy(label_rr,"RESpace");
	sprintf(filename,"%s/Mu_decomp/MD_%c%c_%s_logk_z%05.2f_dim%04d_size%04.0fMpc.dat",MMRRM_PS_OP,label_Aprox_O,label_Aprox_T,label_rr,zrl,dim_nbody,box_size);
}


