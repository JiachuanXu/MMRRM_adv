#include "header/read_args.h"
#include "header/read_write_box.h"
#include "header/cal_Tb.h"

/*
	USAGE: Edit the argumets in "../prmt/mmrrm_adv.h" and parameters in "/header/mmrrm_adv.h", then makefile. Run the executable file ./MMRRM_adv
	WARNING: The output will cover the file with the same name, so make sure to modify the existing file name if needed.
*/

	char boxes_PATH[M_PATH];
	struct realization * info; 
	float box_size;
	int dim_nbody, dim_rt, nr;
	FILE *LOG;

void getname_Tb(int idir, char *filename,double *ave,int real_ct);
void getname_List(char *filename, int nr);
void getname_GloEvo(char *filename, struct tm *now);

int main(int argc, char *argv[])
{
	//Variable declaration
	time_t start_time, curr_time;
	char LOG_NAME[M_BOXNAME],cmd[M_CMD],LIST_NAME[M_BOXNAME];
	int idir;
	struct tm *local;
	FILE *oplist, *GLOB_EVOL;

	start_time=time(NULL);
	local=localtime(&start_time);
	system("mkdir ../Log_files");
	system("mkdir ../Log_files/MMRRM_adv");
	sprintf(LOG_NAME,"../Log_files/MMRRM_adv/MMRRM_adv_log_file_%d_%d_%d_%d:%d:%d",local->tm_year+YEAR_START,local->tm_mon+MON_START,local->tm_mday,local->tm_hour,local->tm_min,local->tm_sec);
	sprintf(cmd,"mkdir %s",MMRRM_BOX_OP);
	system(cmd);
	printf("Creating log file...\n");
	LOG=fopen(LOG_NAME,"a");
	if(LOG==NULL){
		printf("*** MMRRM_adv.c: Error: Can't create log file!\n");
		exit(EXIT_FAILURE);
	}
	else
		fprintf(LOG," Log files for MMRRM_adv\n Date: %s\n",asctime(local));

	printf("Reading arguments...\n");
	fprintf(LOG, "Reading arguments...\n");
	info = read_args(boxes_PATH, &box_size, &dim_nbody,&dim_rt,&nr);
	if(STATIS==1){
		system("mkdir ../Log_files/Output_List");
		/*if(OPTHIN){
			if(HIGHTS)
				sprintf(LIST_NAME,"../Log_files/Output_List/Delta_Tb_List_TH_size%04.0f_dim%04d_count%03d",box_size,dim_nbody,nr);
			else
				sprintf(LIST_NAME,"../Log_files/Output_List/Delta_Tb_List_TN_size%04.0f_dim%04d_count%03d",box_size,dim_nbody,nr);
		}
		else{
			if(HIGHTS)
				sprintf(LIST_NAME,"../Log_files/Output_List/Delta_Tb_List_FH_size%04.0f_dim%04d_count%03d",box_size,dim_nbody,nr);
			else
				sprintf(LIST_NAME,"../Log_files/Output_List/Delta_Tb_List_FN_size%04.0f_dim%04d_count%03d",box_size,dim_nbody,nr);
		}*/
		getname_List(LIST_NAME, nr);
		oplist=fopen(LIST_NAME,"w");
		if(oplist==NULL){
			printf("*** MMRRM_adv.c: Error: Can't create output list file!\n");
			exit(EXIT_FAILURE);
		}
		else
			fprintf(oplist,"%d\n%d\n%f\n%d\n%d\n%d\n",nr,dim_nbody,box_size,OPTHIN,HIGHTS,MESH2MESH);
	}
	if(GLOBAL_EVOL){
		char GE_Name[M_BOXNAME];
		system("mkdir ../Log_files/Global_Evolution");
		/*if(OPTHIN){
			if(HIGHTS)
				sprintf(GE_Name,"../Log_files/Global_Evolution/TH_dim%04d_size%04.0f_%d_%d_%d_%d:%d:%d", dim_nbody, box_size,local->tm_year+YEAR_START,local->tm_mon+MON_START,local->tm_mday,local->tm_hour,local->tm_min,local->tm_sec);
			else
				sprintf(GE_Name,"../Log_files/Global_Evolution/TN_dim%04d_size%04.0f_%d_%d_%d_%d:%d:%d", dim_nbody, box_size,local->tm_year+YEAR_START,local->tm_mon+MON_START,local->tm_mday,local->tm_hour,local->tm_min,local->tm_sec);
		}
		else{
			if(HIGHTS)
				sprintf(GE_Name,"../Log_files/Global_Evolution/FH_dim%04d_size%04.0f_%d_%d_%d_%d:%d:%d", dim_nbody, box_size,local->tm_year+YEAR_START,local->tm_mon+MON_START,local->tm_mday,local->tm_hour,local->tm_min,local->tm_sec);
			else
				sprintf(GE_Name,"../Log_files/Global_Evolution/FN_dim%04d_size%04.0f_%d_%d_%d_%d:%d:%d", dim_nbody, box_size,local->tm_year+YEAR_START,local->tm_mon+MON_START,local->tm_mday,local->tm_hour,local->tm_min,local->tm_sec);
		}*/
		getname_GloEvo(GE_Name, local);
		GLOB_EVOL = fopen(GE_Name,"w");
	}
	float Ts_Evol[nr],Tcmb_Evol[nr],xHIv_Evol[nr],xHIm_Evol[nr],Op_Thick_Evol[nr];
	double Tb_Evol[nr];
	/******************************* Looping through zrl array *******************************/
	for(int real_ct=0;real_ct<nr;real_ct++){
		float * delta;	//Density fluctuation
		float * xHI;	//Neutral fraction
		float * Ts; 	//Spin Temp
		double * Tb;		//Brightness Temp
		float * v;      //Peculiar velocity, cMpc/s
		char op_name[3][M_BOXNAME];
		double Tb_ave[4];
		for(int p=0;p<4;p++){
			Tb_ave[p]=0.0;
		}
		printf("\n****** Processing redshift: %.2f******\n\nReading input files...\n",info[real_ct].zrl);
		fprintf(LOG, "\n****** Processing redshift: %.2f******\n\nReading input files...\n",info[real_ct].zrl);
	
		delta = fmyalloc(box_vol(dim_nbody));
		xHI = fmyalloc(box_vol(dim_rt));
		Ts = fmyalloc(box_vol(dim_nbody));
		fbox_read(delta,boxes_PATH,info[real_ct].n_den,dim_nbody);// Overdensity
		xHIv_Evol[real_ct] = fbox_read(xHI,boxes_PATH,info[real_ct].n_xHI,dim_rt);//neutral fraction
		Ts_Evol[real_ct] = fbox_read(Ts,boxes_PATH,info[real_ct].n_Ts,dim_nbody);
		xHIm_Evol[real_ct] = Mass_Ave(delta, xHI, dim_nbody, dim_rt);
		Tcmb_Evol[real_ct] = T_cmb(info[real_ct].zrl);
		Op_Thick_Evol[real_ct] = 0.0;
		
		/********************* Looping through LoS direction: Calculate Tb *********************/
	
		for(idir=1;idir<4;idir++){
			printf("Looping along %d direction...\n",idir);
			fprintf(LOG, "Looping along %d direction...\n",idir);
			v = fmyalloc(box_vol(dim_nbody));
			Tb = dmyalloc(box_vol(dim_nbody));
			switch(idir){
				case 1:
					fbox_read(v,boxes_PATH,info[real_ct].n_vx,dim_nbody);
					break;
				case 2:
					fbox_read(v,boxes_PATH,info[real_ct].n_vy,dim_nbody);
					break;
				case 3:
					fbox_read(v,boxes_PATH,info[real_ct].n_vz,dim_nbody);
					break;
				default:
					printf("*** MMRRM_adv.c: Error: idir out of range!\n");
					fprintf(LOG,"*** MMRRM_adv.c: Error: idir out of range!\n");
					fclose(LOG);
					exit(EXIT_FAILURE);
			}
			
			Op_Thick_Evol[real_ct] += cal_Tb(delta,xHI,Ts,v,Tb,idir,info[real_ct].zrl,Tb_ave + idir);
			free(v); v = NULL;
			printf("           Tb average: %e\n",Tb_ave[idir]);
			fprintf(LOG,"           Tb average: %e\n",Tb_ave[idir]);
			printf("Writing the Tb%d box...\n",idir);
			fprintf(LOG, "Writing the Tb%d box...\n",idir);
			getname_Tb(idir, op_name[idir-1],Tb_ave,real_ct);
			dbox_write(Tb,op_name[idir-1],dim_nbody);
			free(Tb); Tb = NULL;
		}	
		/********************* 				LoS direction: End			 *********************/	
		Tb_ave[0]=(Tb_ave[1]+Tb_ave[2]+Tb_ave[3])/3.0;
		if(GLOBAL_EVOL){
			Tb_Evol[real_ct] = Tb_ave[0];
			Op_Thick_Evol[real_ct] /= 3.0;
		}
		free(delta); delta = NULL;
		free(xHI); xHI = NULL;
		free(Ts); Ts = NULL;
		if(POWINK||POWINLOGK){
			sprintf(cmd,"./test_powers %s %s %s",op_name[0],op_name[1],op_name[2]);
			system(cmd);
		}
		if(TBPDF){
			sprintf(cmd,"./test_pdf %s %s %s",op_name[0],op_name[1],op_name[2]);
			system(cmd);
		}
		if(MUDECINK||MUDECINLOGK){
			sprintf(cmd,"./test_mudecomp %s %s %s %s %s %s %s %d",op_name[0],op_name[1],op_name[2], boxes_PATH, info[real_ct].n_den, info[real_ct].n_xHI, info[real_ct].n_Ts, dim_rt);
			system(cmd);
		}
		
		if(STATIS==1)
			fprintf(oplist,"%05.2f\n%s\n%s\n%s\n",info[real_ct].zrl,op_name[0],op_name[1],op_name[2]);
		curr_time=time(NULL);
		printf("\tTime elapsed: %.2f minutes\n",difftime(curr_time,start_time)/60.0);
		fprintf(LOG, "\tTime elapsed: %.2f minutes\n",difftime(curr_time,start_time)/60.0);
	}
	/******************************* zrl array: End *******************************/
	if(STATIS==1){
		fclose(oplist);
		sprintf(cmd,"./test_statool %s",LIST_NAME);
		system(cmd);
	}
	if(GLOBAL_EVOL){
		for(int i=0;i<nr;i++)
			fprintf(GLOB_EVOL, "%.2f\t%e\t%e\t%e\t%e\t%e\t%f\n", info[i].zrl, Ts_Evol[i],Tcmb_Evol[i],Tb_Evol[i],xHIv_Evol[i],xHIm_Evol[i],Op_Thick_Evol[i]); // redshift, Spin Temperature, T_CMB, 21 cm Brightness Temperature, <xHI>_cell, <xHI>_mass, Optical Thick Cells' Fraction
		fclose(GLOB_EVOL);
	}
	free(info); info = NULL;
	curr_time=time(NULL);
	printf("\tTotal time elapsed: %.2f minutes\n",difftime(curr_time,start_time)/60.0);
	fprintf(LOG, "\tTotal time elapsed: %.2f minutes\n",difftime(curr_time,start_time)/60.0);
	if(fclose(LOG)!=0){
		printf("*** MMRRM_adv.c: Error: Can't close log file!\n");
		exit(EXIT_FAILURE);
	}
	return 0;
}

void getname_Tb(int idir, char *filename,double *ave, int real_ct){
	char label_idir, label_Aprox_O, label_Aprox_T;
	char label_rr[8] = {'0'};
	(HIGHTS)? (label_Aprox_T = 'H'):(label_Aprox_T = 'N');
	(OPTHIN)? (label_Aprox_O = 'T'):(label_Aprox_O = 'F');
	(MESH2MESH)? strcpy(label_rr, "RSSpace"):strcpy(label_rr,"RESpace");
	switch(idir){
		case 1:
			label_idir = 'x';
			break;
		case 2:
			label_idir = 'y';
			break;
		case 3:
			label_idir = 'z';
			break;
		default:
		printf("*** MMRRM_adv.c: Error: idir out of range!\n");
		fprintf(LOG,"*** MMRRM_adv.c: Error: idir out of range!\n");
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	sprintf(filename,"Tb%c_Aprox%c%c_%s_z%05.2f_size%04.0fMpc_dim%04d_Ave%e",label_idir, label_Aprox_O, label_Aprox_T, label_rr, info[real_ct].zrl, box_size, dim_nbody, *((double *)ave + idir));
	/*switch(idir){
		case 1:
			if(OPTHIN==0){//Op-thin limit
				if(HIGHTS)
					sprintf(op_name,"Tbx_FH_z%05.2f_size%04.0fMpc_dim%04d_Ave%e",info[real_ct].zrl,box_size,dim_nbody,*((double *)ave + idir));
				else
					sprintf(op_name,"Tbx_FN_z%05.2f_size%04.0fMpc_dim%04d_Ave%e",info[real_ct].zrl,box_size,dim_nbody,*((double *)ave + idir));
			}
			else{
				if(HIGHTS)
					sprintf(op_name,"Tbx_TH_z%05.2f_size%04.0fMpc_dim%04d_Ave%e",info[real_ct].zrl,box_size,dim_nbody,*((double *)ave + idir));
				else
					sprintf(op_name,"Tbx_TN_z%05.2f_size%04.0fMpc_dim%04d_Ave%e",info[real_ct].zrl,box_size,dim_nbody,*((double *)ave + idir));
			}
			break;
		case 2:
			if(OPTHIN==0){
				if(HIGHTS)
					sprintf(op_name,"Tby_FH_z%05.2f_size%04.0fMpc_dim%04d_Ave%e",info[real_ct].zrl,box_size,dim_nbody,*((double *)ave + idir));
				else
					sprintf(op_name,"Tby_FN_z%05.2f_size%04.0fMpc_dim%04d_Ave%e",info[real_ct].zrl,box_size,dim_nbody,*((double *)ave + idir));
			}
			else{
				if(HIGHTS)
					sprintf(op_name,"Tby_TH_z%05.2f_size%04.0fMpc_dim%04d_Ave%e",info[real_ct].zrl,box_size,dim_nbody,*((double *)ave + idir));
				else
					sprintf(op_name,"Tby_TN_z%05.2f_size%04.0fMpc_dim%04d_Ave%e",info[real_ct].zrl,box_size,dim_nbody,*((double *)ave + idir));
			}
			break;
		case 3:
			if(OPTHIN==0){
				if(HIGHTS)
					sprintf(op_name,"Tbz_FH_z%05.2f_size%04.0fMpc_dim%04d_Ave%e",info[real_ct].zrl,box_size,dim_nbody,*((double *)ave + idir));
				else
					sprintf(op_name,"Tbz_FN_z%05.2f_size%04.0fMpc_dim%04d_Ave%e",info[real_ct].zrl,box_size,dim_nbody,*((double *)ave + idir));
			}
			else{
				if(HIGHTS)
					sprintf(op_name,"Tbz_TH_z%05.2f_size%04.0fMpc_dim%04d_Ave%e",info[real_ct].zrl,box_size,dim_nbody,*((double *)ave + idir));
				else
					sprintf(op_name,"Tbz_TN_z%05.2f_size%04.0fMpc_dim%04d_Ave%e",info[real_ct].zrl,box_size,dim_nbody,*((double *)ave + idir));
			}
			break;
		default:
			printf("Error: idir out of range!\n");
			fprintf(LOG,"Error: idir out of range!\n");
			fclose(LOG);
			exit(EXIT_FAILURE);
	}*/
}

void getname_List(char *filename, int nr){
	char label_Aprox_O, label_Aprox_T;
	char label_rr[8] = {'\0'};
	(HIGHTS)? (label_Aprox_T = 'H'):(label_Aprox_T = 'N');
	(OPTHIN)? (label_Aprox_O = 'T'):(label_Aprox_O = 'F');
	(MESH2MESH)? strcpy(label_rr, "RSSpace"):strcpy(label_rr,"RESpace"); // RedShifted Space & REal Space
	sprintf(filename,"../Log_files/Output_List/Delta_Tb_List_%c%c_%s_size%04.0f_dim%04d_count%03d", label_Aprox_O, label_Aprox_T, label_rr,box_size,dim_nbody,nr);
}

void getname_GloEvo(char *filename, struct tm *now){
	char label_Aprox_O, label_Aprox_T;
	char label_rr[8] = {'\0'};
	(HIGHTS)? (label_Aprox_T = 'H'):(label_Aprox_T = 'N');
	(OPTHIN)? (label_Aprox_O = 'T'):(label_Aprox_O = 'F');
	(MESH2MESH)? strcpy(label_rr, "RSSpace"):strcpy(label_rr,"RESpace");
	sprintf(filename,"../Log_files/Global_Evolution/%c%c_%s_dim%04d_size%04.0f_%d_%d_%d_%d:%d:%d", label_Aprox_O, label_Aprox_T, label_rr, dim_nbody, box_size,now->tm_year+YEAR_START,now->tm_mon+MON_START,now->tm_mday,now->tm_hour,now->tm_min,now->tm_sec);
}

