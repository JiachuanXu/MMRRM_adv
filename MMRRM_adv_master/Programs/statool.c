#include "header/read_args.h"
#include "header/read_write_box.h"
#define STATBIN 50 // Number of bins of skewness and kurtosis
/*
	Usage: ./statool <PATH: List name>
	Output format: redshift, mean, deviation, rms, skewness, kurtosis
*/

	FILE *LOG;
	int nr, dim, opthin, hights, mesh2mesh;
	float box_size;

void getname_stat(char *filename, struct tm *now);

int main(int argc, char *argv[]){
	FILE *LIST,*STAT;
	char c_nr[M_FLOAT],c_z[M_FLOAT],c_dim[M_FLOAT],c_size[M_FLOAT],c_flagO[M_FLOAT],c_flagT[M_FLOAT],c_flagmm[M_FLOAT],STAT_NAME[M_BOXNAME],LOG_NAME[M_BOXNAME],cmd[M_CMD];
	struct TbList *list;
	time_t start_time;
	struct tm *local;

	start_time=time(NULL);
	local=localtime(&start_time);
	system("mkdir ../Log_files");
	system("mkdir ../Log_files/Statistics");
	/*if(OPTHIN){
		if(HIGHTS)
			sprintf(STAT_NAME,"%s/Statistics/Stat_info_TH_%d_%d_%d_%d:%d:%d",MMRRM_PS_OP,local->tm_year+YEAR_START,local->tm_mon+MON_START,local->tm_mday,local->tm_hour,local->tm_min,local->tm_sec);
		else
			sprintf(STAT_NAME,"%s/Statistics/Stat_info_TN_%d_%d_%d_%d:%d:%d",MMRRM_PS_OP,local->tm_year+YEAR_START,local->tm_mon+MON_START,local->tm_mday,local->tm_hour,local->tm_min,local->tm_sec);
	}
	else{
		if(HIGHTS)
			sprintf(STAT_NAME,"%s/Statistics/Stat_info_FH_%d_%d_%d_%d:%d:%d",MMRRM_PS_OP,local->tm_year+YEAR_START,local->tm_mon+MON_START,local->tm_mday,local->tm_hour,local->tm_min,local->tm_sec);
		else
			sprintf(STAT_NAME,"%s/Statistics/Stat_info_FN_%d_%d_%d_%d:%d:%d",MMRRM_PS_OP,local->tm_year+YEAR_START,local->tm_mon+MON_START,local->tm_mday,local->tm_hour,local->tm_min,local->tm_sec);
	}*/
	sprintf(LOG_NAME,"../Log_files/Statistics/Statool_log_%d_%d_%d_%d:%d:%d",local->tm_year+YEAR_START,local->tm_mon+MON_START,local->tm_mday,local->tm_hour,local->tm_min,local->tm_sec);
	LOG=fopen(LOG_NAME,"w");

	if(argc!=2){
		printf("USAGE: ./test_statool <Box list>\n");
		exit(EXIT_FAILURE);
	}
	printf("****** Calling statistics tools");
	printf("Reading arguments...\n");
	fprintf(LOG,"Reading arguments...\n");
	LIST=fopen(argv[1],"r");
	arg_fgets(c_nr,M_FLOAT,LIST);
	arg_fgets(c_dim,M_FLOAT,LIST);
	arg_fgets(c_size,M_FLOAT,LIST);
	arg_fgets(c_flagT,M_FLOAT,LIST);
	arg_fgets(c_flagO,M_FLOAT,LIST);
	arg_fgets(c_flagmm,M_FLOAT,LIST);
	nr = atoi(c_nr);
	dim = atoi(c_dim);
	box_size = atof(c_size);
	opthin = atoi(c_flagO);
	hights = atoi(c_flagT);
	mesh2mesh = atoi(c_flagmm);

	list=(struct TbList *)malloc(sizeof(struct TbList)*nr);
	for(int ct=0;ct<nr;ct++) {
		arg_fgets(c_z,M_FLOAT,LIST);
		(*((struct TbList *)list+ct)).zrl = atof(c_z);
		arg_fgets((*((struct TbList *)list+ct)).fTb_x,M_BOXNAME,LIST);
		arg_fgets((*((struct TbList *)list+ct)).fTb_y,M_BOXNAME,LIST);
		arg_fgets((*((struct TbList *)list+ct)).fTb_z,M_BOXNAME,LIST);
	}
	fclose(LIST);
	sprintf(cmd,"mkdir %s",MMRRM_PS_OP);
	system(cmd);
	sprintf(cmd,"mkdir %s/Statistics",MMRRM_PS_OP);
	system(cmd);
	getname_stat(STAT_NAME, local);
	STAT=fopen(STAT_NAME,"w");
	fprintf(STAT,"%5s\t%12s\t%12s\t%12s\t%12s\t%12s\n","zrl","Tb_mean","Tb_dev","Tb_rms","Skewness","Kurtosis");
	
/************************* Start Calculating Statistics ************************/

	for(int real_ct=0;real_ct<nr;real_ct++){
		printf("Looping through realization %d...\n",real_ct+1);
		fprintf(LOG,"Looping through realization %d...\n",real_ct+1);
		double mean[4],rms[4],skew[4],kurt[4],dev[4];
		double *my_mean,*my_rms,*my_skew,*my_kurt,*my_dev;
		for(int i=0;i<4;i++){
			mean[i] = 0.0;
			rms[i] = 0.0;
			skew[i] = 0.0;
			kurt[i] = 0.0;
			dev[i] = 0.0;
		}
		for(int idir=0;idir<3;idir++){
			double *Tb;
			Tb=(double *)malloc(sizeof(double)*box_vol(dim));
			switch(idir){
				case 0:
				dbox_read(Tb,MMRRM_BOX_OP,list[real_ct].fTb_x,dim);
				break;
				case 1:
				dbox_read(Tb,MMRRM_BOX_OP,list[real_ct].fTb_y,dim);
				break;
				case 2:
				dbox_read(Tb,MMRRM_BOX_OP,list[real_ct].fTb_z,dim);
				break;
				default:
				exit(EXIT_FAILURE);
			}
			unsigned long long ct;
			my_mean = dmyalloc(nthreads);
			my_rms = dmyalloc(nthreads);
			my_skew = dmyalloc(nthreads);
			my_kurt = dmyalloc(nthreads);
			my_dev = dmyalloc(nthreads);

			#pragma omp parallel num_threads(nthreads) default(none) shared(idir,dim,my_mean,my_rms,mean,rms,Tb) private(ct)
			{
				int my_rank = omp_get_thread_num();
				#pragma omp for
				for(ct=0;ct<box_vol(dim);ct++){
					my_mean[my_rank] += Tb[ct];
					my_rms[my_rank] += Tb[ct] * Tb[ct]; 
				}
				#pragma omp atomic
				mean[idir] += my_mean[my_rank];
				#pragma omp atomic
				rms[idir] += my_rms[my_rank];
			}
			mean[idir] /= (double)box_vol(dim);
			rms[idir] /= (double)box_vol(dim);
			#pragma omp parallel num_threads(nthreads) default(none) shared(idir,dim,my_dev,my_skew,my_kurt,dev,skew,kurt,Tb, mean) private(ct)
			{
				int my_rank = omp_get_thread_num();
				#pragma omp for
				for(ct=0;ct<box_vol(dim);ct++){
					double momentum = (Tb[ct] - mean[idir]);
					my_dev[my_rank] += pow(momentum,2.0);
					my_skew[my_rank] += pow(momentum,3.0);
					my_kurt[my_rank] += pow(momentum,4.0);
				}
				#pragma omp atomic
				dev[idir] += my_dev[my_rank];
				#pragma omp atomic
				skew[idir] += my_skew[my_rank];
			 	#pragma omp atomic
				kurt[idir] += my_kurt[my_rank];
			}
			dev[idir] /= (double)(box_vol(dim)-1.0);
			dev[idir] = sqrt(dev[idir]);
			skew[idir] /= (double)box_vol(dim)*pow(dev[idir],3.0);
			kurt[idir] /= (double)box_vol(dim)*pow(dev[idir],4.0);
			free(Tb); Tb = NULL;
		}
		mean[0] = (mean[1]+mean[2]+mean[3])/3.0;
		dev[0] = (dev[1]+dev[2]+dev[3])/3.0;
		rms[0] = (rms[1]+rms[2]+rms[3])/3.0;
		skew[0] = (skew[1]+skew[2]+skew[3])/3.0;
		kurt[0] = (kurt[1]+kurt[2]+kurt[3])/3.0;
		fprintf(STAT,"%5.2f\t%e\t%e\t%e\t%e\t%e\n",list[real_ct].zrl,mean[0],dev[0],rms[0],skew[0],kurt[0]);
	}
	fclose(STAT);
	free(list); list = NULL;
	printf("Normal ending...\n");
	fprintf(LOG,"Normal ending...\n");
	fclose(LOG);

	return 0;
}

void getname_stat(char *filename, struct tm *now){
	char label_Aprox_O, label_Aprox_T;
	char label_rr[8] = {'\0'};
	(hights)? (label_Aprox_T = 'H'):(label_Aprox_T = 'N');
	(opthin)? (label_Aprox_O = 'T'):(label_Aprox_O = 'F');
	(mesh2mesh)? strcpy(label_rr, "RSSpace"):strcpy(label_rr,"RESpace");
	sprintf(filename,"%s/Statistics/Stat_info_%c%c_%s_dim%d_size%4.0fMpc_%d_%d_%d_%d:%d:%d",MMRRM_PS_OP,label_Aprox_O,label_Aprox_T,label_rr,dim,box_size,now->tm_year+YEAR_START,now->tm_mon+MON_START,now->tm_mday,now->tm_hour,now->tm_min,now->tm_sec);
}
