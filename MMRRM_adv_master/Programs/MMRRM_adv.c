// MMRRM_adv.c 
// Main program 
// Generate $\delta T_b$ datacube utilizing 
// mass overdensity, peculiar velocity, spin temperature, neutral fraction
// Conduct statictical analysis
// Edit history:
// 2018.4.20: Add more comments, Jiachuan Xu

/*
	USAGE: Edit the argumets in 
	        "../prmt/mmrrm_adv.h" 
	       and parameters in 
	        "/header/mmrrm_adv.h" 
	       Then makefile. Run the executable file ./MMRRM_adv
	WARNING: The output will cover the file with the same name.
*/

#include "header/read_args.h"
#include "header/read_write_box.h"
#include "header/cal_Tb.h"

/*******************	Variable Declaration	*******************/

char boxes_PATH[M_PATH];
struct realization * info;// contents for data cubes
float box_size;
int dim_nbody, dim_rt;// dimension of N-body & RT
int nr;// number of data cubes to be processed
FILE *LOG;

/********************	Function Declaration	********************/

// getname_Tb: generate the filename for $\delta T_b$ data cube
// Inputs:
//		idir:	LoS direction, 1-x 2-y 3-z
//		ave:	average $\delta Tb$, [mK]
//		real_ct:id for the data cube
// Output:
// 		filename
void getname_Tb(int idir, char *filename,double *ave,int real_ct);

// getname_List: generate filename for $\delta Tb$ boxes list 
// Inputs:
//		nr:		number of data cubes
// Outputs:
//		filename
void getname_List(char *filename, int nr);

// getname_GloEvo: generate filename for global evolution log file 
// Inputs:
//		tm:		[structure] time
// Outputs:
//		filename
void getname_GloEvo(char *filename, struct tm *now);

/***********************	Main Program	***********************/

int main(int argc, char *argv[])
{
	//Variable declaration
	time_t start_time, curr_time;
	char LOG_NAME[M_BOXNAME],cmd[M_CMD],LIST_NAME[M_BOXNAME];
	int idir;
	struct tm *local;
	FILE *oplist, *GLOB_EVOL;
	
	/********************** Initialization **********************/
	// Time initialization
	start_time=time(NULL);
	local=localtime(&start_time);
	// Create folders
	system("mkdir ../Log_files");
	system("mkdir ../Log_files/MMRRM_adv");
	// Creat log file 
	sprintf(LOG_NAME,"../Log_files/MMRRM_adv/MMRRM_adv_log_file_\
%d_%d_%d_%d:%d:%d",local->tm_year+YEAR_START,local->tm_mon+MON_START,
local->tm_mday,local->tm_hour,local->tm_min,local->tm_sec);// Name of log file <----Log file name
	sprintf(cmd,"mkdir %s",MMRRM_BOX_OP);
	system(cmd);
	printf("Creating log file...\n");
	LOG=fopen(LOG_NAME,"a");
	if(LOG==NULL){
		printf("*** MMRRM_adv.c: Error: Can't create log file!\n");
		exit(EXIT_FAILURE);
	}
	else
		fprintf(LOG," Log files for MMRRM_adv\n Date: %s\n"
			,asctime(local));
	
	// Read arguments in ../prmt
	printf("Reading arguments...\n");
	fprintf(LOG, "Reading arguments...\n");
	// the information of data cube, like component fields' filename,
	// redshift, are stored in structure "realization"
	info = read_args(boxes_PATH, &box_size, &dim_nbody,&dim_rt,&nr);
	// If want to calculate kurtosis, skewness, <rms>...
	// create list recording output boxes' filename
	if(STATIS==1){
		system("mkdir ../Log_files/Output_List");
		getname_List(LIST_NAME, nr);
		oplist=fopen(LIST_NAME,"w");
		if(oplist==NULL){
			printf("*** MMRRM_adv.c: \
Error: Can't create output list file!\n");
			exit(EXIT_FAILURE);
		}
		else// dump basic info 
			fprintf(oplist,"%d\n%d\n%f\n%d\n%d\n%d\n",
				nr,dim_nbody,box_size,OPTHIN,HIGHTS,MESH2MESH);
	}
	// If want to calculate global evolution of T_CMB, Ts, $\delta T_b$
	// <xHI>_v, <xHI>_m and optical thick cells 
	// Create global evolution log file 
	if(GLOBAL_EVOL){
		char GE_Name[M_BOXNAME];
		system("mkdir ../Log_files/Global_Evolution");
		getname_GloEvo(GE_Name, local);
		GLOB_EVOL = fopen(GE_Name,"w");
	}

	// global evolution variables
	float Ts_Evol[nr],Tcmb_Evol[nr]; 
	float xHIv_Evol[nr],xHIm_Evol[nr];
	float Op_Thick_Evol[nr];// calculate fraction, so float
	double Tb_Evol[nr];

	/************************** Contents **************************/

	// Start for: data cubes (realizations)
	for(int real_ct=0;real_ct<nr;real_ct++)
	{
		float * delta;	// mass overdensity [dimensionless]
		float * xHI;	// neutral fraction [dimensionless]
		float * Ts; 	// spin temperature [K]
		double * Tb;	// brightness temperature difference [mK or K]
		float * v;      // peculiar velocity, [cMpc/s]
		char op_name[3][M_BOXNAME]; // output $\delta T_b$ filename
		double Tb_ave[4];// <$\delta T_b$>
		for(int p=0;p<4;p++){
			Tb_ave[p]=0.0;
		}
		Op_Thick_Evol[real_ct] = 0.0; // optical-thick cell counts
		
		// Read-in component fields
		printf("\n****** Processing redshift: %.2f******\n\n\
Reading input files...\n",info[real_ct].zrl);
		fprintf(LOG, "\n****** Processing redshift: %.2f******\n\n\
Reading input files...\n",info[real_ct].zrl);
		// memory allocation (peculiar velocity will be allocated later)
		delta = fmyalloc(box_vol(dim_nbody));
		xHI = fmyalloc(box_vol(dim_rt));
		Ts = fmyalloc(box_vol(dim_nbody));
		// read-in & record global info of <xHI>_v, <xHI>_m, Ts, T_CMB
		fbox_read(delta,boxes_PATH,info[real_ct].n_den,dim_nbody);
		xHIv_Evol[real_ct] = 
		fbox_read(xHI,boxes_PATH,info[real_ct].n_xHI,dim_rt);
		xHIm_Evol[real_ct] = 
		Mass_Ave(delta, xHI, dim_nbody, dim_rt);
		// If force Ts >> T_CMB and optical-thin (HIGHTS=1 OPTHIN==1)
		// do not read-in Ts data cube and feed program with fake Ts <--------------fake Ts flag
		if (HIGHTS==1 && OPTHIN==1){
			Ts_Evol[real_ct] = FLT_MAX;
			for (unsigned long long i; i<box_vol(dim_nbody), i++){
				Ts[i] = FLT_MAX;
			}
		}
		else{
			Ts_Evol[real_ct] = 
			fbox_read(Ts,boxes_PATH,info[real_ct].n_Ts,dim_nbody);
		}
		Tcmb_Evol[real_ct] = 
		T_cmb(info[real_ct].zrl);
		
		/********* Looping through LoS direction: Calculate Tb *********/
	
		// Start for: LoS
		for(idir=1;idir<4;idir++)
		{
			printf("Looping along %d direction...\n",idir);
			fprintf(LOG, "Looping along %d direction...\n",idir);
			
			// allocate memory for LoS peculiar velocity & output
			v = fmyalloc(box_vol(dim_nbody));
			Tb = dmyalloc(box_vol(dim_nbody));
			// read-in 
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
					fprintf(LOG,
						"*** MMRRM_adv.c: Error: idir out of range!\n");
					fclose(LOG);
					exit(EXIT_FAILURE);
			}
			
			// calculate $\delta T_b$
			//          &
			// record number of optical-thick cells
			Op_Thick_Evol[real_ct] += 
			cal_Tb(delta,xHI,Ts,v,Tb,idir,info[real_ct].zrl,Tb_ave + idir);

			// get output filename & dump $\delta T_b$ into that
			free(v); v = NULL;
			printf("           Tb average: %e\n",Tb_ave[idir]);
			fprintf(LOG,"           Tb average: %e\n",Tb_ave[idir]);
			printf("Writing the Tb%d box...\n",idir);
			fprintf(LOG, "Writing the Tb%d box...\n",idir); 
			getname_Tb(idir, op_name[idir-1],Tb_ave,real_ct);
			dbox_write(Tb,op_name[idir-1],dim_nbody);
			free(Tb); Tb = NULL;
		}// End for: LoS

		// calculate LoS-averaged <$\delta T_b$>
		Tb_ave[0]=(Tb_ave[1]+Tb_ave[2]+Tb_ave[3])/3.0;
		// write into global info 
		if(GLOBAL_EVOL){
			Tb_Evol[real_ct] = Tb_ave[0];
			Op_Thick_Evol[real_ct] /= 3.0;
		}
		free(delta); delta = NULL;
		free(xHI); xHI = NULL;
		free(Ts); Ts = NULL;
		// If want to calculate power spectrum of $\delta T_b$
		// call subroutine test_powers
		if(POWINK||POWINLOGK){
			sprintf(cmd,"./test_powers %s %s %s",
				op_name[0],op_name[1],op_name[2]);
			system(cmd);
		}
		// If want to calculate $\delta T_b$ PDF 
		// call subroutine test_pdf
		if(TBPDF){
			sprintf(cmd,"./test_pdf %s %s %s",
				op_name[0],op_name[1],op_name[2]);
			system(cmd);
		}
		// If want to calculate mu-decomposition
		// call subroutine test_mudecomp 
		if(MUDECINK||MUDECINLOGK){
			sprintf(cmd,"./test_mudecomp %s %s %s %s %s %s %s %d",
				op_name[0],op_name[1],op_name[2], boxes_PATH, 
				info[real_ct].n_den, info[real_ct].n_xHI, 
				info[real_ct].n_Ts, dim_rt);
			system(cmd);
		}
		// If want to calculate kurtosis, skewness, <rms>...
		// write $\delta T_b$ filename infos
		if(STATIS==1)
			fprintf(oplist,"%05.2f\n%s\n%s\n%s\n",
				info[real_ct].zrl,op_name[0],op_name[1],op_name[2]);

		// timing
		curr_time=time(NULL);
		printf("\tTime elapsed: %.2f minutes\n",
			difftime(curr_time,start_time)/60.0);
		fprintf(LOG, "\tTime elapsed: %.2f minutes\n",
			difftime(curr_time,start_time)/60.0);
	}// End for: data cubes (realizations)

	// If want to calculate kurtosis, skewness, <rms>...
	// call subroutine test_statool
	if(STATIS==1){
		fclose(oplist);
		sprintf(cmd,"./test_statool %s",LIST_NAME);
		system(cmd);
	}
	// If want to calculate global evolution
	// write global info:
	// z, Ts, T_CMB, $\delta T_b$, <xHI>_v, <xHI>_m, f_op-thick
	// 1,  K,     K,      mK or K,       1,       1,          1
	if(GLOBAL_EVOL){
		for(int i=0;i<nr;i++)
			fprintf(GLOB_EVOL, "%.2f\t%e\t%e\t%e\t%e\t%e\t%f\n", 
			info[i].zrl, Ts_Evol[i],Tcmb_Evol[i],Tb_Evol[i],
			xHIv_Evol[i],xHIm_Evol[i],Op_Thick_Evol[i]); 
		fclose(GLOB_EVOL);
	}
	// clear up
	free(info); info = NULL;
	curr_time=time(NULL);
	printf("\tTotal time elapsed: %.2f minutes\n",
		difftime(curr_time,start_time)/60.0);
	fprintf(LOG, "\tTotal time elapsed: %.2f minutes\n",
		difftime(curr_time,start_time)/60.0);
	if(fclose(LOG)!=0){
		printf("*** MMRRM_adv.c: Error: Can't close log file!\n");
		exit(EXIT_FAILURE);
	}
	return 0;
}

void getname_Tb(int idir, char *filename,double *ave, int real_ct)
{
	// label for approximations in the filename:
	// force optical-thin & force high Ts: TH 
	// force optical-thin & ~force high Ts: TN 
	// ~force optical-thin & force high Ts: FH 
	// ~force optical-thin & ~force high Ts: FN [benchmark]
	// $\delta T_b$ observed in real space: RE
	// $\delta T_b$ observed in redshift space: RS 
	// LoS: x or y or z 
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
	sprintf(filename,
		"Tb%c_Aprox%c%c_%s_z%05.2f_size%04.0fMpc_dim%04d_Ave%e",
		label_idir, label_Aprox_O, label_Aprox_T, label_rr, 
		info[real_ct].zrl, box_size, dim_nbody, 
		*((double *)ave + idir));
}

void getname_List(char *filename, int nr)
{
	// label for approximations in the filename:
	// force optical-thin & force high Ts: TH 
	// force optical-thin & ~force high Ts: TN 
	// ~force optical-thin & force high Ts: FH 
	// ~force optical-thin & ~force high Ts: FN [benchmark]
	// $\delta T_b$ observed in real space: RE
	// $\delta T_b$ observed in redshift space: RS 
	// count: number of data cubes
	char label_Aprox_O, label_Aprox_T;
	char label_rr[8] = {'\0'};
	(HIGHTS)? (label_Aprox_T = 'H'):(label_Aprox_T = 'N');
	(OPTHIN)? (label_Aprox_O = 'T'):(label_Aprox_O = 'F');
	(MESH2MESH)? strcpy(label_rr, "RSSpace"):strcpy(label_rr,"RESpace");
	sprintf(filename,"../Log_files/Output_List/\
Delta_Tb_List_%c%c_%s_size%04.0f_dim%04d_count%03d", 
label_Aprox_O, label_Aprox_T, label_rr,box_size,dim_nbody,nr);
}

void getname_GloEvo(char *filename, struct tm *now)
{
	char label_Aprox_O, label_Aprox_T;
	char label_rr[8] = {'\0'};
	(HIGHTS)? (label_Aprox_T = 'H'):(label_Aprox_T = 'N');
	(OPTHIN)? (label_Aprox_O = 'T'):(label_Aprox_O = 'F');
	(MESH2MESH)? strcpy(label_rr, "RSSpace"):strcpy(label_rr,"RESpace");
	sprintf(filename, "../Log_files/Global_Evolution/\
%c%c_%s_dim%04d_size%04.0f_%d_%d_%d_%d:%d:%d", 
label_Aprox_O, label_Aprox_T, label_rr, dim_nbody, box_size,
now->tm_year+YEAR_START,now->tm_mon+MON_START,now->tm_mday,
now->tm_hour,now->tm_min,now->tm_sec);
}