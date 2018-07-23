// pixelfluxpdf.c
// calculate the probability distribution function of
// delta Tb
// Edit history:
// 2018.7.16: Add more comments

/*
	USAGE: ./pdf <Delta Tb file x> <Delta Tb file y> <Delta Tb file z> [pdf file]
	Output format: difference probability
*/

#include "header/read_write_box.h"
#include "header/read_args.h"
#define PDFBIN 100

/******************************************************************/
/*******************	Variable Declaration	*******************/
/******************************************************************/

	FILE *LOG;
	float zrl,box_size;
	int dim, mesh2mesh, opthin, hights;

/******************************************************************/
/********************	Function Declaration	*******************/
/******************************************************************/

// getname_PDF: generate the filename of the pdf file
// Output:
// 		filename
void getname_PDF(char *filename);

/******************************************************************/
/***********************	Main Program	***********************/
/******************************************************************/

int main(int argc, char *argv[])
{

	// Variable declaration
	FILE *PDF;
	double *Tb;
	double count_bin[PDFBIN],Tb_diff_bin[PDFBIN];
	char PDF_NAME[M_BOXNAME],LOG_NAME[M_BOXNAME],cmd[M_CMD];
	time_t start_time;
	struct tm *local;
	start_time=time(NULL);
	local=localtime(&start_time);

	/********************** Initialization **********************/
	// arguments check
	if(argc!=4&&argc!=5){
		printf("USAGE: ./test_pdf <Delta Tb file x> <Delta Tb file y> \
<Delta Tb file z> [pdf file]\n");
		exit(EXIT_FAILURE);
	}
	// open log file
	if(access(MMRRM_PS_OP,F_OK)==-1){
		sprintf(cmd,"mkdir %s",MMRRM_PS_OP);
		system(cmd);
	}
	char MT[M_BOXNAME];
	sprintf(MT,"%s/Pixelflux_PDF",MMRRM_PS_OP);
	if(access(MT,F_OK)==-1){
		sprintf(cmd,"mkdir %s/Pixelflux_PDF",MMRRM_PS_OP);
		system(cmd);
	}
	if(access("../Log_files",F_OK)==-1)
		system("mkdir ../Log_files");
	if(access("../Log_files/PDF",F_OK)==-1)
		system("mkdir ../Log_files/PDF");
	sprintf(LOG_NAME,"../Log_files/PDF/PDF_log_%d_%d_%d_%d:%d:%d",
		local->tm_year+YEAR_START,local->tm_mon+MON_START,local->tm_mday,
		local->tm_hour,local->tm_min,local->tm_sec);
	LOG=fopen(LOG_NAME,"w");

	// read parameters from filename
	read_title(argv[1],&zrl,&box_size,&dim,&opthin,&hights,&mesh2mesh);
	printf("****** Calling PDF: z = %.2f\n",zrl);
	for(int ct=0;ct<PDFBIN;ct++){
		count_bin[ct] = 0.0;
		Tb_diff_bin[ct] = 0.0;
	}
	double ave[4],max,min,dif;
	Tb = (double *)malloc(sizeof(double)*box_vol(dim)*3);
	// Initializing bins
	// set the mean, max and the min of delta Tb
	for (int idir = 1;idir < 4; idir ++){
		dbox_read(Tb + box_vol(dim)*(idir - 1),MMRRM_BOX_OP,argv[idir],dim);

		ave[idir] = 0.0;
		if(idir==1){
			max = Tb[0];
			min = Tb[0];
		}
		for(unsigned long long ct=0;ct<box_vol(dim);ct++){
			if(Tb[ct + (idir - 1)*box_vol(dim)]>max)
				max = Tb[ct + (idir - 1)*box_vol(dim)];
			if(Tb[ct + (idir -1)*box_vol(dim)]<min)
				min = Tb[ct + (idir - 1)*box_vol(dim)];
			ave[idir] += Tb[ct + (idir - 1)*box_vol(dim)];
		}
		ave[idir] /= (double)box_vol(dim);
	}
	ave[0] = (ave[1]+ave[2]+ave[3])/3.0;
	dif = (max-min)/(PDFBIN-1);
	// Looping
	for(int idir=1;idir<4;idir++){
		printf("Looping through dir %d\n",idir);
		int i,j,k;
		#pragma omp parallel num_threads(nthreads) default(none) shared(idir,dim,min,dif,count_bin,Tb,Tb_diff_bin) private(i,j,k)
		{
			#pragma omp for
			for(i=0;i<dim;i++){
				for(j=0;j<dim;j++){
					for(k=0;k<dim;k++){
						int ref;
						ref = (int)((Tb[in_tr(i,j,k,dim) + (idir - 1)*box_vol(dim)] - min)/dif + 0.5);
					# pragma omp atomic
						count_bin[ref]++;
					# pragma omp atomic
						Tb_diff_bin[ref] += Tb[in_tr(i,j,k,dim) + (idir - 1)*box_vol(dim)];
					}
				}
			}
		}
	}
	free(Tb); Tb = NULL;
	for(int i = 0; i < PDFBIN; i++){
		if(count_bin[i]!=0)
			Tb_diff_bin[i] /= count_bin[i];
		else
			Tb_diff_bin[i] = min + dif * i;
	}
	// write pdf into file
	if(argc==5)
		sprintf(PDF_NAME,"%s/Pixelflux_PDF/%s",MMRRM_PS_OP,argv[4]);
	else
		getname_PDF(PDF_NAME);
	if((PDF=fopen(PDF_NAME,"w"))==NULL){
		printf("*** pixelfluxpdf.c: Error: Can't creat pdf file\n");
		exit(EXIT_FAILURE);
	}
	// print header
	fprintf(PDF, "# delta_Tb counts_fraction\n");
	for(int ct=0;ct<PDFBIN;ct++){
		count_bin[ct] /= (double)(3.0*box_vol(dim));
		//Tb_diff_bin[ct] -= ave[0];
		fprintf(PDF,"%e\t%e\n",Tb_diff_bin[ct],count_bin[ct]);
	}
	fclose(PDF);
	fclose(LOG);
	printf("Normal Ending...\n");
	return 0;
}

void getname_PDF(char *filename){
	char label_Aprox_O, label_Aprox_T;
	char label_rr[8] = {'\0'};
	(hights)? (label_Aprox_T = 'H'):(label_Aprox_T = 'N');
	(opthin)? (label_Aprox_O = 'T'):(label_Aprox_O = 'F');
	(mesh2mesh)? strcpy(label_rr, "RSSpace"):strcpy(label_rr,"RESpace");
	sprintf(filename,"%s/Pixelflux_PDF/PDF_%c%c_%s_z%05.2f_size%04.0f_dim%04d_bin%d",MMRRM_PS_OP,label_Aprox_O,label_Aprox_T,label_rr,zrl,box_size,dim,PDFBIN);
}
