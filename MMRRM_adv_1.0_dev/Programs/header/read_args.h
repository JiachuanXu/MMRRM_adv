#include "mmrrm_adv.h"

/*
USAGE: To read the arguments, args_mmrrm_adv.dat
*/
extern FILE *LOG;
struct realization{
	char n_den[M_BOXNAME];
	char n_xHI[M_BOXNAME];
	char n_Ts[M_BOXNAME];
	char n_vx[M_BOXNAME];
	char n_vy[M_BOXNAME];
	char n_vz[M_BOXNAME];
	float zrl;
};
struct TbList{
	float zrl;
	char fTb_x[M_BOXNAME];
	char fTb_y[M_BOXNAME];
	char fTb_z[M_BOXNAME];
};

struct realization * read_args(char *boxes_PATH,float *box_size, int *dim_nbody, int *dim_rt, int *nr);
void read_title(char *filename,float *zrl,float *box_size,int *dim, int *flag_O, int *flag_T, int *flag_mm);
void arg_fgets(char *list, int linemax, FILE *f_p);

struct realization * read_args(char *boxes_PATH,float *box_size, int *dim_nbody, int *dim_rt, int *nr)
{
	char c_z[M_FLOAT], c_dim_nbody[M_FLOAT], c_dim_rt[M_FLOAT], c_box_size[M_FLOAT],c_nr[M_FLOAT];
	struct realization * pinfo;
	FILE * fp;
	if ((fp=fopen(std_arg,"r"))==NULL)
	{
		printf("Error: Can't open arguments %s!\n",std_arg);
		fprintf(LOG,"Error: Can't open arguments %s!\n",std_arg);
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	arg_fgets(boxes_PATH,M_PATH,fp);
	arg_fgets(c_dim_nbody,M_FLOAT,fp);
	arg_fgets(c_dim_rt,M_FLOAT,fp);
	arg_fgets(c_box_size,M_FLOAT,fp);
	arg_fgets(c_nr,M_FLOAT,fp);
	*box_size=atof(c_box_size);
	*dim_nbody=atoi(c_dim_nbody);
	*dim_rt = atoi(c_dim_rt);
	*nr=atoi(c_nr);
	pinfo = (struct realization *)malloc(sizeof(struct realization)*(*nr));
	printf("%33s: %s\n","Dir of boxes", boxes_PATH);
	printf("%33s: %.0f\n","Box size", *box_size);
	printf("%33s: %d\n","Dimension of box (nbody)", *dim_nbody);
	printf("%33s: %d\n","Dimension of box (radiation transfer)", *dim_rt);
	printf("%33s: %d\n","Number of realization", *nr);
	for(int ct=0;ct<(*nr);ct++) {
		arg_fgets((*((struct realization *)pinfo+ct)).n_den,M_BOXNAME,fp);
		arg_fgets((*((struct realization *)pinfo+ct)).n_xHI,M_BOXNAME,fp);
		arg_fgets((*((struct realization *)pinfo+ct)).n_Ts,M_BOXNAME,fp);
		arg_fgets((*((struct realization *)pinfo+ct)).n_vx,M_BOXNAME,fp);
		arg_fgets((*((struct realization *)pinfo+ct)).n_vy,M_BOXNAME,fp);
		arg_fgets((*((struct realization *)pinfo+ct)).n_vz,M_BOXNAME,fp);
		arg_fgets(c_z,M_FLOAT,fp);
		(*((struct realization *)pinfo+ct)).zrl = atof(c_z);
	}
	if(fclose(fp)!=0){
		printf("Error: Cant' close file %s\n", std_arg);
		fprintf(LOG,"Error: Cant' close file %s\n", std_arg);
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	if(((*dim_nbody)%(*dim_rt))!=0){
		printf("read_args.h: Error: dim_nbody/dim_rt must be an integer!");
		exit(EXIT_FAILURE);
	}
	return pinfo;
	
}

void read_title(char *filename,float *zrl,float *box_size,int *dim,int *flag_O, int *flag_T, int *flag_mm){
	printf("Reading file title\n");
	fprintf(LOG,"Reading file title\n");
	char *find_z = strstr(filename,"_z");
	char *find_size = strstr(filename,"_size");
	char *find_dim = strstr(filename,"_dim");
	char *find_flagmm = strstr(filename, "Space_");
	char *find_flagOT = strstr(filename,"_Aprox");
	char c_redshift[M_FLOAT] = {'\0'};
	//c_redshift[0]='\0';
	char c_boxsize[M_FLOAT] = {'\0'};
	//c_boxsize[0]='\0';
	char c_dim[M_FLOAT] = {'\0'};
	//c_dim[0]='\0';
	char c_mm[M_FLOAT] = {'\0'};
	char c_flagO[M_FLOAT] = {'\0'};
	char c_flagT[M_FLOAT] = {'\0'};

	if(find_z == NULL){
		printf("Error: Can't read redshift\n");
		exit(1);
	}
	else
		strncat(c_redshift,find_z+2,5);

	if(find_size == NULL){
		printf("Error: Can't read boxsize\n");
		exit(1);
	}
	else
		strncat(c_boxsize,find_size+5,4);

	if(find_dim == NULL){
		printf("Error: Can't read dimension\n");
		exit(1);
	}
	else
		strncat(c_dim,find_dim+4,4);

	if(find_flagmm == NULL){
		printf("Error: Can't read real or redshift space\n");
		exit(1);
	}
	else
		strncat(c_mm,find_flagmm - 2,2);

	if(find_flagOT == NULL){
		printf("Error: Can't read approximation\n");
		exit(1);
	}
	else
		strncat(c_flagO,find_flagOT + 6,1);

	if(find_flagOT == NULL){
		printf("Error: Can't read approximation\n");
		exit(1);
	}
	else
		strncat(c_flagT,find_flagOT + 7,1);


	*zrl = atof(c_redshift);
	*box_size = atof(c_boxsize);
	*dim = atoi(c_dim);

	if(!strcmp(c_mm,"RE"))
		*flag_mm = 0;
	else if (!strcmp(c_mm, "RS"))
		*flag_mm = 1;
	else{
		printf("Error: Can't \"read real or redshift space\"\n");
		exit(1);
	}

	if(!strcmp(c_flagO, "F"))
		*flag_O = 0;
	else if(!strcmp(c_flagO, "T"))
		*flag_O = 1;
	else{
		printf("Error: Can't read \"optical thin or full-treatment\"\n");
		exit(1);
	}

	if(!strcmp(c_flagT, "N"))
		*flag_T = 0;
	else if(!strcmp(c_flagT,"H"))
		*flag_T = 1;
	else{
		printf("Error: Can't read \"High Ts limit or Normal Ts\"\n");
		exit(1);
	}
}

void arg_fgets(char *list, int linemax, FILE *f_p){
	char list_temp[M_TEMP];
	char *find;
	fgets(list_temp,linemax,f_p);
	while(list_temp[0]=='\n')
		fgets(list_temp,linemax,f_p);
	find=strchr(list_temp,'\n');
	if(find==NULL){
		printf("Error: Error happens while reading argument! Maybe the argument lenth run out of the limit defined in mmrrm_adv.h, or you didn't insert line break between two terms.\n");
		fprintf(LOG,"Error: Error happens while reading argument! Maybe the argument lenth run out of the limit defined in mmrrm_adv.h, or you didn't insert line break between two terms.\n");
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	else
		*find='\0';
	strcpy(list,list_temp);
}