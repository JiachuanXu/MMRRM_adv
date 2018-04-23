#include "header/read_write_box.h"
#include "header/read_args.h"
#include "header/cal_fftw.h"
/*
	Usage: ./powers <Tb x> <Tb y> <Tb z>
	Output format:
	linear: kint k_magnitude(1/Mpc) power(mK^2) deviation
	log:    k_magnitude(1/Mpc) power(mK^2) deviation
*/

	FILE *LOG;
	float z,box_size;
	int dim, mesh2mesh, opthin, hights;

void getname_PSLin(char *filename);
void getname_PSLog(char *filename);

int main(int argc, char *argv[]){
	//read the fft box
	char PS_LIN_NAME[M_BOXNAME],PS_LOG_NAME[M_BOXNAME],LOG_NAME[M_BOXNAME],cmd[M_CMD];
	time_t start_time;
	struct tm *local;
	FILE *PS_LIN,*PS_LOG;
	int mid,nintvl;
	double *power_spec_lin,*err_lin;
	double *power_spec_log,*err_log;
	double *k_array_lin, *k_array_log;
	double k_step,dlnk;
	

	if(argc!=4){
		printf("Usage: ./powers <Tb along x> <Tb along y> <Tb along z>\n");
		exit(EXIT_FAILURE);
	}
	start_time=time(NULL);
	local=localtime(&start_time);
	system("mkdir ../Log_files");
	system("mkdir ../Log_files/Powers");
	sprintf(LOG_NAME,"../Log_files/Powers/SAPS_logfile_%d_%d_%d_%d:%d:%d",local->tm_year+YEAR_START,local->tm_mon+MON_START,local->tm_mday,local->tm_hour,local->tm_min,local->tm_sec);
	LOG=fopen(LOG_NAME,"a");
	if(LOG==NULL){
		printf("*** powers.c: Error: Can't create log file!\n");
		exit(EXIT_FAILURE);
	}
	else
		fprintf(LOG," Log files for powers_k\n Date: %s\n",asctime(local));
	sprintf(cmd,"mkdir %s",MMRRM_PS_OP);
	system(cmd);

// Reading file title 
	read_title(argv[1],&z,&box_size,&dim,&opthin,&hights,&mesh2mesh);
	mid=dim/2;
	k_step=2*PI/box_size;
	printf("****** Calling Power Spectrum: z = %.2f\n",z);
	fprintf(LOG,"z:%05.2f\tsize:%04.0f\tdim:%04d\n",z,box_size,dim);
/***************************** Looping through LoS   *******************************/
	int idir;
	double *Tb;
	unsigned long long *num_lin_bin,*num_log_bin;
	fftw_complex *power;

	if(POWINK){
		power_spec_lin = dmyalloc((unsigned long long)(mid-1));
		err_lin = dmyalloc((unsigned long long)(mid-1));
		num_lin_bin = ullmyalloc((unsigned long long)(mid-1));
		k_array_lin = dmyalloc(mid - 1);
	}
	if(POWINLOGK){
		nintvl = (int)(log((double)(mid-1)/(double)kmincut)/log(((double)kmincut+1.0)/(double)kmincut));
		dlnk = (double)(log((double)(mid-1)/(double)kmincut)/((double)nintvl-1.0));
		power_spec_log = dmyalloc((unsigned long long)nintvl);
		err_log = dmyalloc((unsigned long long)nintvl);
		num_log_bin = ullmyalloc((unsigned long long)nintvl);
		k_array_log = dmyalloc((unsigned long long)nintvl);
	}

	for(idir=1;idir<4;idir++){
		Tb = dmyalloc(box_vol(dim));
		dbox_read(Tb,MMRRM_BOX_OP,argv[idir],dim);
		
		if((power = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*kbox_vol(dim)))==NULL){
			printf("*** powers.c: Error: Can't allocate meomory for fft box\n");
			fprintf(LOG,"*** powers.c: Error: Can't allocate memory for fft box\n");
			fclose(LOG);
			exit(1);
		}
		calculate_power(Tb,Tb,dim,power);
		free(Tb); Tb = NULL;
		
		if(POWINK)
			spherical_sorting_linear(power, num_lin_bin, power_spec_lin, err_lin, dim);
		if(POWINLOGK)
			spherical_sorting_log(power, num_log_bin, power_spec_log, err_log, dim, nintvl, dlnk);
		fftw_free(power); power = NULL;
	}
	if(POWINK){
		spherical_average_linear(num_lin_bin, power_spec_lin, err_lin, k_array_lin, dim, mid, k_step);
		free(num_lin_bin); num_lin_bin = NULL;
	}
	if(POWINLOGK){
		spherical_average_log(num_log_bin, power_spec_log, err_log, k_array_log, dim, nintvl, dlnk, k_step);
		free(num_log_bin); num_log_bin = NULL;
	}

	if(POWINK){
		sprintf(cmd,"mkdir %s/SAPS_lin",MMRRM_PS_OP);
		system(cmd);
		/*if(OPTHIN==1){
			if(HIGHTS)
				sprintf(PS_LIN_NAME,"%s/SAPS_lin/SAPS_TH_k_z%05.2f_dim%04d_size%04.0fMpc.dat",MMRRM_PS_OP,z,dim,box_size);
			else
				sprintf(PS_LIN_NAME,"%s/SAPS_lin/SAPS_TN_k_z%05.2f_dim%04d_size%04.0fMpc.dat",MMRRM_PS_OP,z,dim,box_size);
		}
		else{
			if(HIGHTS)
				sprintf(PS_LIN_NAME,"%s/SAPS_lin/SAPS_FH_k_z%05.2f_dim%04d_size%04.0fMpc.dat",MMRRM_PS_OP,z,dim,box_size);
			else
				sprintf(PS_LIN_NAME,"%s/SAPS_lin/SAPS_FN_k_z%05.2f_dim%04d_size%04.0fMpc.dat",MMRRM_PS_OP,z,dim,box_size);
		}*/
		getname_PSLin(PS_LIN_NAME);
		if((PS_LIN=fopen(PS_LIN_NAME,"w"))==NULL){
			printf("*** powers.c: Error: Can't open power spectrum file %s\n",PS_LIN_NAME);
			fprintf(LOG,"*** powers.c: Error: Can't open power spectrum file %s\n",PS_LIN_NAME);
			fclose(LOG);
			exit(EXIT_FAILURE);
		}
		fprintf(PS_LIN, "# k\tps\tps_err\n");
		for(int i=0;i<mid-1;i++)
			fprintf(PS_LIN,"%E\t%E\t%E\n",k_array_lin[i],power_spec_lin[i],err_lin[i]*power_spec_lin[i]);// kint k_magnitude power deviation
		fclose(PS_LIN);
		free(power_spec_lin); power_spec_lin = NULL;
		free(err_lin); err_lin = NULL;
		free(k_array_lin); k_array_lin = NULL;
	}
	if(POWINLOGK){
		sprintf(cmd,"mkdir %s/SAPS_log",MMRRM_PS_OP);
		system(cmd);
		/*if(OPTHIN==1){
			if(HIGHTS)
				sprintf(PS_LOG_NAME,"%s/SAPS_log/SAPS_TH_logk_z%05.2f_dim%04d_size%04.0fMpc.dat",MMRRM_PS_OP,z,dim,box_size);
			else
				sprintf(PS_LOG_NAME,"%s/SAPS_log/SAPS_TN_logk_z%05.2f_dim%04d_size%04.0fMpc.dat",MMRRM_PS_OP,z,dim,box_size);
		}
		else{
			if(HIGHTS)
				sprintf(PS_LOG_NAME,"%s/SAPS_log/SAPS_FH_logk_z%05.2f_dim%04d_size%04.0fMpc.dat",MMRRM_PS_OP,z,dim,box_size);
			else
				sprintf(PS_LOG_NAME,"%s/SAPS_log/SAPS_FN_logk_z%05.2f_dim%04d_size%04.0fMpc.dat",MMRRM_PS_OP,z,dim,box_size);
		}*/
		getname_PSLog(PS_LOG_NAME);
		if((PS_LOG=fopen(PS_LOG_NAME,"w"))==NULL){
			printf("*** powers.c: Error: Can't open power spectrum file %s\n",PS_LOG_NAME);
			fprintf(LOG,"*** powers.c: Error: Can't open power spectrum file %s\n",PS_LOG_NAME);
			fclose(LOG);
			exit(EXIT_FAILURE);
		}
		fprintf(PS_LOG, "# k\tps\tps_err\n");
		for(int q=0;q<nintvl;q++)
			fprintf(PS_LOG,"%E\t%E\t%E\n",k_array_log[q],power_spec_log[q],power_spec_log[q]*err_log[q]);// k_magnitude power deviation
		fclose(PS_LOG);
		free(power_spec_log); power_spec_log = NULL;
		free(err_log); err_log = NULL;
		free(k_array_log); k_array_log = NULL;
	}
	printf("Normal ending...\n");
	fprintf(LOG,"Normal ending...\n");
	fclose(LOG);
	return 0;
}

void getname_PSLin(char *filename){
	char label_Aprox_O, label_Aprox_T;
	char label_rr[8] = {'\0'};
	(hights)? (label_Aprox_T = 'H'):(label_Aprox_T = 'N');
	(opthin)? (label_Aprox_O = 'T'):(label_Aprox_O = 'F');
	(mesh2mesh)? strcpy(label_rr, "RSSpace"):strcpy(label_rr,"RESpace");
	sprintf(filename,"%s/SAPS_lin/SAPS_%c%c_%s_k_z%05.2f_dim%04d_size%04.0fMpc.dat",MMRRM_PS_OP,label_Aprox_O,label_Aprox_T,label_rr,z,dim,box_size);
}
void getname_PSLog(char *filename){
	char label_Aprox_O, label_Aprox_T;
	char label_rr[8] = {'\0'};
	(hights)? (label_Aprox_T = 'H'):(label_Aprox_T = 'N');
	(opthin)? (label_Aprox_O = 'T'):(label_Aprox_O = 'F');
	(mesh2mesh)? strcpy(label_rr, "RSSpace"):strcpy(label_rr,"RESpace");
	sprintf(filename,"%s/SAPS_log/SAPS_%c%c_%s_logk_z%05.2f_dim%04d_size%04.0fMpc.dat",MMRRM_PS_OP,label_Aprox_O,label_Aprox_T,label_rr,z,dim,box_size);
}

