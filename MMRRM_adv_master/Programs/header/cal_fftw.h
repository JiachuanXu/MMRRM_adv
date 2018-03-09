#include "mmrrm_adv.h"
//#include "cal_Tb.h"
#define ZERO_MEAN 0 // Transform into zero-mean box

FILE *LOG;
void calculate_power(double *in_1,double *in_2,int dim,fftw_complex *power);
double fluct(double *in,int dim);
void spherical_sorting_linear(fftw_complex *power, unsigned long long *num_lin_bin, double * power_spec_lin, double * err_lin, int dim);
void spherical_sorting_log(fftw_complex * power,unsigned long long *num_log_bin, double * power_spec_log, double * err_log, int dim, int nintvl, double dlnk);
void spherical_average_linear(unsigned long long *num_lin_bin, double *power_spec_lin, double *err_lin, double *k_array_lin, int dim, int mid, double k_step);
void spherical_average_log(unsigned long long *num_log_bin, double *power_spec_log, double *err_log, double *k_array_log, int dim, int nintvl, double dlnk, double k_step);
void azimuthal_sorting(fftw_complex *power, unsigned long long *count_decomp, double *D2_decomp, double *sqr_decomp, int dim_nbody, int mid, int idir);
void azimuthal_average(unsigned long long *count_decomp, double *D2_decomp, double *sqr_decomp, double *mu_decomp, int mid);


void calculate_power(double *in_1,double *in_2,int dim,fftw_complex *power){
	fftw_plan p;
	fftw_complex *out_1,*out_2;
	printf("Calculating powers...\n");
	if(in_1!=in_2){
		if((out_1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*kbox_vol(dim)))==NULL){
			printf("Error: Can't allocate meomory for fft box\n");
			fprintf(LOG,"Error: Can't allocate memory for fft box\n");
			fclose(LOG);
			exit(EXIT_FAILURE);
		}
		if((out_2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*kbox_vol(dim)))==NULL){
			printf("Error: Can't allocate meomory for fft box\n");
			fprintf(LOG,"Error: Can't allocate memory for fft box\n");
			fclose(LOG);
			exit(EXIT_FAILURE);
		}
		if(ZERO_MEAN){
			fluct(in_1,dim);
			fluct(in_2,dim);
		}
		p=fftw_plan_dft_r2c_3d(dim,dim,dim,in_1,out_1,FFTW_ESTIMATE);
		fftw_execute(p);
		fftw_destroy_plan(p);
		p=fftw_plan_dft_r2c_3d(dim,dim,dim,in_2,out_2,FFTW_ESTIMATE);
		fftw_execute(p);
		fftw_destroy_plan(p);
		for(unsigned long long ct=0;ct<kbox_vol(dim);ct++)
		power[ct] = out_1[ct]*conj(out_2[ct]);
		fftw_free(out_1);
		fftw_free(out_2);
	}
	else{
		if((out_1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*kbox_vol(dim)))==NULL){
			printf("Error: Can't allocate meomory for fft box\n");
			fprintf(LOG,"Error: Can't allocate memory for fft box\n");
			fclose(LOG);
			exit(EXIT_FAILURE);
		}
		if(ZERO_MEAN){
			fluct(in_1,dim);
		}
		p=fftw_plan_dft_r2c_3d(dim,dim,dim,in_1,out_1,FFTW_ESTIMATE);
		fftw_execute(p);
		fftw_destroy_plan(p);
		for(unsigned long long ct=0;ct<kbox_vol(dim);ct++)
			power[ct] = out_1[ct]*conj(out_1[ct]);
		fftw_free(out_1);
	}
}

double fluct(double *in,int dim){
	double ave=0.0,new_ave=0.0;
	for(int ct=0;ct<box_vol(dim);ct++)
		ave += in[ct];
	ave /= pow((double)dim,3.0);
	for(int ct=0;ct<box_vol(dim);ct++){
		in[ct] /= ave;
		in[ct] -= (double)1.0;
		new_ave += in[ct];
	}
	new_ave /= pow((double)dim,3.0);
	printf("Average: %e\nDelta average: %e\n",ave,new_ave);
	fprintf(LOG,"Average: %e\nDelta average: %e\n",ave,new_ave);
	return ave;
}

void spherical_sorting_linear(fftw_complex *power,unsigned long long *num_lin_bin, double *power_spec_lin, double *err_lin, int dim){
	int i,j,k,ki,kj,kk,kint_lin, mid = dim/2;
	#pragma omp parallel num_threads(nthreads) default(none) shared(dim, mid,power,num_lin_bin,power_spec_lin,err_lin) private(i,j,k,ki,kj,kk,kint_lin)
	{
		# pragma omp for
		for(i=0;i<dim;i++){
			if(i<(mid+1))
				ki=i;
			else
				ki=i-dim;
			for(j=0;j<dim;j++){
				if(j<(mid+1))
				kj=j;
				else
					kj=j-dim;
				for(k=0;k<(mid+1);k++){
					kk=k;

					kint_lin=(int)(sqrt(pow((double)ki,2.0)+pow((double)kj,2.0)+pow((double)kk,2.0))+0.5);
					if((kint_lin>0)&&(kint_lin<mid)&&(ki!=mid)&&(kj!=mid)&&(kk!=mid)){
						if(kk==0){
						# pragma omp atomic
							num_lin_bin[kint_lin-1]++;
						# pragma omp atomic
								power_spec_lin[kint_lin-1] += power[in_tr_c(i,j,k,dim)];
						# pragma omp atomic
						err_lin[kint_lin-1] += power[in_tr_c(i,j,k,dim)]*power[in_tr_c(i,j,k,dim)];
						}
						else{
						# pragma omp atomic
							num_lin_bin[kint_lin-1]+=2;
						# pragma omp atomic
							power_spec_lin[kint_lin-1] += power[in_tr_c(i,j,k,dim)]*2;
						# pragma omp atomic
							err_lin[kint_lin-1] += power[in_tr_c(i,j,k,dim)]*power[in_tr_c(i,j,k,dim)]*2;
						}
					}	
				}
			}
		}
	}
}

void spherical_sorting_log(fftw_complex *power,unsigned long long *num_log_bin, double *power_spec_log, double *err_log, int dim, int nintvl, double dlnk){
	int i,j,k,ki,kj,kk,kint_log,mid = dim/2;
	#pragma omp parallel num_threads(nthreads) default(none) shared(dim, mid,power,dlnk,nintvl,num_log_bin,power_spec_log,err_log) private(i,j,k,ki,kj,kk,kint_log)
	{
		# pragma omp for
		for(i=0;i<dim;i++){
			if(i<(mid+1))
				ki=i;
			else
				ki=i-dim;
			for(j=0;j<dim;j++){
				if(j<(mid+1))
				kj=j;
				else
					kj=j-dim;
				for(k=0;k<(mid+1);k++){
					kk=k;
				
					kint_log=(int)(log(sqrt(pow((double)ki,2.0)+pow((double)kj,2.0)+pow((double)kk,2.0))/(double)kmincut)/dlnk + 0.5);
					if((kint_log>=0)&&(kint_log<nintvl)&&(ki!=mid)&&(kj!=mid)&&(kk!=mid)){
						if(kk==0){
						# pragma omp atomic
							num_log_bin[kint_log]++;
						# pragma omp atomic
							power_spec_log[kint_log] += power[in_tr_c(i,j,k,dim)];
						# pragma omp atomic
							err_log[kint_log] += power[in_tr_c(i,j,k,dim)]*power[in_tr_c(i,j,k,dim)];
						}
						else{
						# pragma omp atomic
							num_log_bin[kint_log]+=2;
						# pragma omp atomic
							power_spec_log[kint_log] += power[in_tr_c(i,j,k,dim)]*2;
						# pragma omp atomic
							err_log[kint_log] += power[in_tr_c(i,j,k,dim)]*power[in_tr_c(i,j,k,dim)]*2;
						}
					}
				}
			}
		}
	}
}

void spherical_average_linear(unsigned long long *num_lin_bin, double *power_spec_lin, double *err_lin, double *k_array_lin, int dim, int mid, double k_step){
	for(int q=0;q<mid-1;q++){
		power_spec_lin[q] /= (double)num_lin_bin[q];
		err_lin[q] /= (double)num_lin_bin[q];
		err_lin[q] = sqrt(err_lin[q]/pow(power_spec_lin[q],2.0)-1.0);
		power_spec_lin[q] = (power_spec_lin[q]*4.0*PI*pow((double)(q+1.0),3.0)/pow((double)dim,6.0));
		k_array_lin[q] = (q + 1.0)*k_step;
	}
}
void spherical_average_log(unsigned long long *num_log_bin, double *power_spec_log, double *err_log, double *k_array_log, int dim, int nintvl, double dlnk, double k_step){
	for(int q=0;q<nintvl;q++){
		power_spec_log[q] /= (double)num_log_bin[q];
		err_log[q] /= (double)num_log_bin[q];
		err_log[q] = sqrt(err_log[q]/pow(power_spec_log[q],2.0)-1.0);
		k_array_log[q] = exp(dlnk*(double)q+log((double)kmincut));
		power_spec_log[q] = (power_spec_log[q]*4.0*PI*pow(k_array_log[q],3.0)/pow((double)dim,6.0));
		k_array_log[q] *= k_step;
	}
}

void azimuthal_sorting(fftw_complex *power, unsigned long long *count_decomp, double *D2_decomp, double *sqr_decomp, int dim_nbody, int mid, int idir){
	int i,j,k,ki,kj,kk,kint_lin,imu;

	#pragma omp parallel num_threads(nthreads) default(none) shared(dim_nbody,mid,idir,power,count_decomp,D2_decomp,sqr_decomp) private(i,j,k,ki,kj,kk,kint_lin,imu)
	{
		#pragma omp for
		for(i=0;i<dim_nbody;i++){
			if(i<(mid+1))
				ki=i;
			else
				ki=i-dim_nbody;
			for(j=0;j<dim_nbody;j++){
				if(j<(mid+1))
				kj=j;
				else
					kj=j-dim_nbody;
				for(k=0;k<(mid+1);k++){
					kk=k;
					switch(idir){
						case 1:
						imu = abs(ki);
						break;
						case 2:
						imu = abs(kj);
						break;
						case 3:
						imu = abs(kk);
						break;
						default:
						printf("Sincerely appreciation for my darling, Gao Shang~\nThank you for you espouse!");
						exit(EXIT_FAILURE);
					}
					kint_lin=(int)(sqrt(pow((double)ki,2.0)+pow((double)kj,2.0)+pow((double)kk,2.0))+0.5);
					if((kint_lin>=(ORDER/2))&&(kint_lin<mid)&&(ki!=mid)&&(kj!=mid)&&(kk!=mid)){
						if(kk==0){
						# pragma omp atomic
							count_decomp[bin_index(kint_lin,imu)]++;
						# pragma omp atomic
						D2_decomp[bin_index(kint_lin,imu)] += power[in_tr_c(i,j,k,dim_nbody)];
						# pragma omp atomic
								sqr_decomp[bin_index(kint_lin,imu)] += power[in_tr_c(i,j,k,dim_nbody)]*power[in_tr_c(i,j,k,dim_nbody)];
						}
					else{
						# pragma omp atomic
							count_decomp[bin_index(kint_lin,imu)] += 2;
						# pragma omp atomic
							D2_decomp[bin_index(kint_lin,imu)] += power[in_tr_c(i,j,k,dim_nbody)]*2;
						# pragma omp atomic
							sqr_decomp[bin_index(kint_lin,imu)] += power[in_tr_c(i,j,k,dim_nbody)]*power[in_tr_c(i,j,k,dim_nbody)]*2;
						}
					}
				}
			}
		}
	}
}
void azimuthal_average(unsigned long long *count_decomp, double *D2_decomp, double *sqr_decomp, double *mu_decomp, int mid){
	for(int l = ORDER/2;l < mid;l++){
		for(int m = 0; m <= l;m++){
			D2_decomp[bin_index(l,m)] /= count_decomp[bin_index(l,m)]; // <D21>
			sqr_decomp[bin_index(l,m)] /= count_decomp[bin_index(l,m)];
			sqr_decomp[bin_index(l,m)] = sqrt(sqr_decomp[bin_index(l,m)] - D2_decomp[bin_index(l,m)]*D2_decomp[bin_index(l,m)]); // deviation \sigma (if you need to change the definition of \sigma!)
			mu_decomp[bin_index(l,m)] = (double)m/(double)l; // mu
		}
	}
}