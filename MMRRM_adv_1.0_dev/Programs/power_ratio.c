#include "header/read_write_box.h"
#include "header/read_args.h"
#include "header/cal_fftw.h"

/*
	USAGE: ./power_ratio <Tb1 along x> <Tb1 along y> <Tb1 along z> <Tb2 along x> <Tb2 along y> <Tb2 along z> [label]
	Output: k[1/Mpc] PS1/PS2 error
*/
	FILE *LOG;
	float z,box_size;
	int dim, mesh2mesh, opthin, hights;

void getname_Ratio(char *filename, int flag, char **cont);

int main(int argc,char *argv[]){

	FILE *fp;
	char RATIO[M_BOXNAME], LOG_NAME[M_BOXNAME], cmd[M_CMD];
	int mid, nintvl, idir;
	//int kint_log;
	double k_step, dlnk;
	//double *power_spec_lin_1, *err_lin_1, *k_array_lin_1;
	//double *power_spec_lin_2, *err_lin_2, *k_array_lin_2;
	//double *ratio_power_spec, *err_compound, *ratio_power_spec_log, *err_compound_log, *k_array;
	double *power_spec_log_1, *err_log_1, *k_array_log_1;
	double *power_spec_log_2, *err_log_2, *k_array_log_2;
	double *ratio_power_spec_log, *err_ratio;
	double *Tb1, *Tb2;
	fftw_complex *power1, *power2;
	unsigned long long *num_log_bin1, *num_log_bin2;
	time_t start_time;
	struct tm *local;
	if((argc != 7)&&(argc != 8)){
		printf("Usage: ./power_ratio <Tb1 along x> <Tb1 along y> <Tb1 along z> <Tb2 along x> <Tb2 along y> <Tb2 along z> [label]\n");
		exit(EXIT_FAILURE);
	}
	start_time=time(NULL);
	local=localtime(&start_time);
	system("mkdir ../Log_files");
	system("mkdir ../Log_files/Power_Ratio");
	sprintf(LOG_NAME,"../Log_files/Power_Ratio/log_%d_%d_%d_%d:%d:%d",local->tm_year+YEAR_START,local->tm_mon+MON_START,local->tm_mday,local->tm_hour,local->tm_min,local->tm_sec);
	LOG=fopen(LOG_NAME,"a");
	if(LOG==NULL){
		printf("*** power_ratio.c: Error: Can't create log file!\n");
		exit(EXIT_FAILURE);
	}
	else
		fprintf(LOG," Log files for power_ratio\n Date: %s\n",asctime(local));
	sprintf(cmd,"mkdir %s",MMRRM_PS_OP);
	system(cmd);
	read_title(argv[1],&z,&box_size,&dim,&opthin,&hights,&mesh2mesh);
	printf("****** Calling power_ratio z = %.2f\n", z);
	mid=dim/2;
	k_step=2*PI/box_size;

	nintvl = (int)(log((double)(mid-1)/(double)kmincut)/log(((double)kmincut+1.0)/(double)kmincut));
	dlnk = (double)(log((double)(mid-1)/(double)kmincut)/((double)nintvl-1.0));
	
	power_spec_log_1 = dmyalloc(nintvl);
	power_spec_log_2 = dmyalloc(nintvl);
	ratio_power_spec_log = dmyalloc(nintvl);
	err_log_1 = dmyalloc(nintvl);
	err_log_2 = dmyalloc(nintvl);
	err_ratio = dmyalloc(nintvl);
	k_array_log_1 = dmyalloc(nintvl);
	k_array_log_2 = dmyalloc(nintvl);
	num_log_bin1 = ullmyalloc(nintvl);
	num_log_bin2 = ullmyalloc(nintvl);

	for(idir=1;idir<4;idir++){
		Tb1 = dmyalloc(box_vol(dim));
		Tb2 = dmyalloc(box_vol(dim));
		dbox_read(Tb1,MMRRM_BOX_OP,argv[idir],dim);
		dbox_read(Tb2,MMRRM_BOX_OP,argv[3 + idir],dim);
		
		if((power1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*kbox_vol(dim)))==NULL){
			printf("*** power_ratio.c: Error: Can't allocate meomory for fft box\n");
			fprintf(LOG,"*** power_ratio.c: Error: Can't allocate memory for fft box\n");
			fclose(LOG);
			exit(1);
		}
		if((power2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*kbox_vol(dim)))==NULL){
			printf("*** power_ratio.c: Error: Can't allocate meomory for fft box\n");
			fprintf(LOG,"*** power_ratio.c: Error: Can't allocate memory for fft box\n");
			fclose(LOG);
			exit(1);
		}
		calculate_power(Tb1,Tb1,dim,power1);
		calculate_power(Tb2,Tb2,dim,power2);
		free(Tb1); Tb1 = NULL;
		free(Tb2); Tb2 = NULL;
		
		spherical_sorting_log(power1, num_log_bin1, power_spec_log_1, err_log_1, dim, nintvl, dlnk);
		spherical_sorting_log(power2, num_log_bin2, power_spec_log_2, err_log_2, dim, nintvl, dlnk);
		fftw_free(power1); power1 = NULL;
		fftw_free(power2); power2 = NULL;
	}
	spherical_average_log(num_log_bin1, power_spec_log_1, err_log_1, k_array_log_1, dim, nintvl, dlnk, k_step);
	spherical_average_log(num_log_bin2, power_spec_log_2, err_log_2, k_array_log_2, dim, nintvl, dlnk, k_step);
	free(num_log_bin1); num_log_bin1 = NULL;
	free(num_log_bin2); num_log_bin2 = NULL;
	
	for(int i = 0; i < nintvl; i++){
		ratio_power_spec_log[i] = power_spec_log_1[i] / power_spec_log_2[i];
		err_ratio[i] = sqrt(err_log_1[i]*err_log_1[i] + err_log_2[i]*err_log_2[i]);
	}
	free(power_spec_log_1); power_spec_log_1 = NULL;
	free(power_spec_log_2); power_spec_log_2 = NULL;
	free(err_log_1); err_log_1 = NULL;
	free(err_log_2); err_log_2 = NULL;


	/*if(argc == 7)
		sprintf(RATIO,"%s/PS_Ratio/Ratio_z%05.2f_dim%04d_size%04.0fMpc", MMRRM_PS_OP, z, dim, box_size);
	else
		sprintf(RATIO,"%s/PS_Ratio/Ratio_z%05.2f_dim%04d_size%04.0fMpc_%s", MMRRM_PS_OP, z, dim, box_size,argv[7]);*/
	sprintf(cmd, "mkdir %s/PS_Ratio", MMRRM_PS_OP);
	system(cmd);
	getname_Ratio(RATIO, argc, argv);
	if((fp = fopen(RATIO, "w"))==NULL){
		printf("*** power_ratio.c: ERROR: Can't creat file\n");
		fprintf(LOG, "*** power_ratio.c: ERROR: Can't creat file\n");
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i < nintvl; i++)
		fprintf(fp,"%e\t%e\t%e\n",k_array_log_1[i], ratio_power_spec_log[i], ratio_power_spec_log[i] * err_ratio[i]);
	fclose(fp);
	free(ratio_power_spec_log);ratio_power_spec_log = NULL;
	free(err_ratio); err_ratio = NULL;
	free(k_array_log_1); k_array_log_1 = NULL;
	free(k_array_log_2); k_array_log_2 = NULL;

	printf("Normal Ending...\n");
	fprintf(LOG,"Normal Ending...\n");
	fclose(LOG);
	return 0;
}

void getname_Ratio(char *filename, int flag, char **cont){
	char label_rr[8] = {'\0'};
	(mesh2mesh)? strcpy(label_rr, "RSSpace"):strcpy(label_rr,"RESpace");
	if(flag == 8)
		sprintf(filename,"%s/PS_Ratio/Ratio_%s_z%05.2f_dim%04d_size%04.0fMpc_%s", MMRRM_PS_OP, label_rr, z, dim, box_size, cont[7]);
	else
		sprintf(filename,"%s/PS_Ratio/Ratio_%s_z%05.2f_dim%04d_size%04.0fMpc", MMRRM_PS_OP, label_rr, z, dim, box_size);
}
