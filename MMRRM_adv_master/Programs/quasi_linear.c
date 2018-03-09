#include "header/cal_Tb.h"
#include "header/cal_fftw.h"
#include "header/read_write_box.h"
#include "header/read_args.h"

float zrl, box_size;
int dim_nbody, dim_rt, mid_nbody;
FILE *LOG;


int main(int argc, char *argv[]){
	int i,j,k,ki,kj,kk,kint_log,nintvl;
	float *Ts, *delta, *xHI;
	double k_step, dlnk;
	double *d_delta, *power_spec, *power_scat, *k_array;
	//double *P_mu4;
	fftw_complex * power;
	unsigned long long *num_log_bin;
	char LOG_NAME[M_BOXNAME], cmd[M_CMD], boxes_PATH[M_PATH], PS_LOG_NAME[M_BOXNAME];
	FILE * PS_LOG;
	struct realization * info;
	int nr;
	time_t start_time;
	struct tm *local;
	if(argc!=4){
		printf("Usage: ./quasi-linear <delta> <xHI> <Ts>\n");
		exit(EXIT_FAILURE);
	}
	start_time=time(NULL);
	local=localtime(&start_time);
	system("mkdir ../Log_files");
	system("mkdir ../Log_files/Quasi_Linear");
	sprintf(LOG_NAME,"../Log_files/Quasi_Linear/log_%d_%d_%d_%d:%d:%d",local->tm_year+YEAR_START,local->tm_mon+MON_START,local->tm_mday,local->tm_hour,local->tm_min,local->tm_sec);
	LOG=fopen(LOG_NAME,"a");
	if(LOG==NULL){
		printf("Error: Can't create log file!\n");
		exit(EXIT_FAILURE);
	}
	else
		fprintf(LOG," Log files for quasi-linear mu decomposition\n Date: %s\n",asctime(local));
	sprintf(cmd,"mkdir %s",MMRRM_PS_OP);
	system(cmd);

	info = read_args(boxes_PATH, &box_size, &dim_nbody,&dim_rt,&nr);
	mid_nbody = dim_nbody/2;
	nintvl = (int)(log((double)(mid_nbody-1)/(double)kmincut)/log(((double)kmincut+1.0)/(double)kmincut));
	dlnk = (double)(log((double)(mid_nbody-1)/(double)kmincut)/((double)nintvl-1.0));
	k_step = 2 * PI / box_size;

	for(int real_ct=0; real_ct< nr; real_ct++){
		double Tbhmean, etamean;
		float xHImmean;
		delta = fmyalloc(box_vol(dim_nbody));
		xHI = fmyalloc(box_vol(dim_rt));
		Ts = fmyalloc(box_vol(dim_nbody));
		d_delta = dmyalloc(box_vol(dim_nbody));

		if((power = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*kbox_vol(dim_nbody)))==NULL){
			printf("Error: Can't allocate meomory for fft box\n");
			fprintf(LOG,"Error: Can't allocate memory for fft box\n");
			fclose(LOG);
			exit(1);
		}

		fbox_read(delta,boxes_PATH,info[real_ct].n_den,dim_nbody);// Overdensity
		fbox_read(xHI,boxes_PATH,info[real_ct].n_xHI,dim_rt);//neutral fraction
		fbox_read(Ts,boxes_PATH,info[real_ct].n_Ts,dim_nbody);
		for(unsigned long long ct = 0; ct < box_vol(dim_nbody); ct++)
			d_delta[ct] = (double)delta[ct];

		calculate_power(d_delta,d_delta,dim_nbody,power);
		xHImmean = Mass_Ave(delta,xHI,dim_nbody,dim_rt);
		Tbhmean = Tb_hat(info[real_ct].zrl,xHImmean);
		etamean = Eta_Ave(Ts, info[real_ct].zrl, dim_nbody);
		free(xHI); xHI = NULL;
		free(Ts); Ts = NULL;
		free(delta); delta = NULL;
		free(d_delta); d_delta = NULL;

		power_spec = dmyalloc(nintvl);
		power_scat = dmyalloc(nintvl);
		num_log_bin = ullmyalloc(nintvl);

		#pragma omp parallel num_threads(nthreads) default(none) shared(dim_nbody,mid_nbody,power,dlnk,nintvl,num_log_bin,power_spec,power_scat) private(i,j,k,ki,kj,kk,kint_log)
		{
			# pragma omp for
			for(i=0;i<dim_nbody;i++){
				if(i<(mid_nbody+1))
					ki = i;
				else
					ki = i - dim_nbody;
				for(j=0;j<dim_nbody;j++){
					if(j<(mid_nbody+1))
						kj=j;
					else
						kj= j - dim_nbody;
					for(k=0;k<(mid_nbody+1);k++){
						kk = k;
						kint_log=(int)(log(sqrt(pow((double)ki,2.0)+pow((double)kj,2.0)+pow((double)kk,2.0))/(double)kmincut)/dlnk + 0.5);
						if((kint_log>=0)&&(kint_log<nintvl)&&(ki!=mid_nbody)&&(kj!=mid_nbody)&&(kk!=mid_nbody)){
							if(kk==0){
							# pragma omp atomic
								num_log_bin[kint_log]++;
							# pragma omp atomic
								power_spec[kint_log] += power[in_tr_c(i,j,k,dim_nbody)];
							# pragma omp atomic
								power_scat[kint_log] += power[in_tr_c(i,j,k,dim_nbody)]*power[in_tr_c(i,j,k,dim_nbody)];
							}
							else{
							# pragma omp atomic
								num_log_bin[kint_log]+=2;
							# pragma omp atomic
								power_spec[kint_log] += power[in_tr_c(i,j,k,dim_nbody)]*2;
							# pragma omp atomic
								power_scat[kint_log] += power[in_tr_c(i,j,k,dim_nbody)]*power[in_tr_c(i,j,k,dim_nbody)]*2;
							}
						}
					}
				}
			}
		}
		free(power); power = NULL;
		for(int q=0;q<nintvl;q++){
			power_spec[q] /= (double)num_log_bin[q];
			power_scat[q] /= (double)num_log_bin[q];
			power_scat[q] = sqrt(power_scat[q]/pow(power_spec[q],2.0)-1.0);
			k_array[q] = exp(dlnk*(double)q+log((double)kmincut));
			power_spec[q] = (pow((Tbhmean * etamean),2.0)* power_spec[q]*4.0*PI*pow(k_array[q],3.0)/pow((double)dim_nbody,6.0));
			k_array[q] *= k_step;
		}
		free(num_log_bin); num_log_bin = NULL;

		if((PS_LOG=fopen(PS_LOG_NAME,"w"))==NULL){
			printf("Error: Can't open power spectrum file %s\n",PS_LOG_NAME);
			fprintf(LOG,"Error: Can't open power spectrum file %s\n",PS_LOG_NAME);
			fclose(LOG);
			exit(EXIT_FAILURE);
		}
		for(int q=0;q<nintvl;q++)
			fprintf(PS_LOG,"%E\t%E\t%E\n",k_array[q],power_spec[q],power_spec[q] * power_scat[q]);// k_magnitude power deviation
		fclose(PS_LOG);
		free(power_spec); power_spec = NULL;
		free(power_scat); power_scat = NULL;
		free(k_array); k_array = NULL;	
	}
	printf("Normal ending\n");
	fprintf(LOG,"Normal ending\n");
	fclose(LOG);
	return 0;

}