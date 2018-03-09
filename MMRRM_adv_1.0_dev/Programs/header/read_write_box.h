#include "mmrrm_adv.h"

/*
USAGE: To read the simulation boxes, including density field, neutral fraction, Ts and peculiar velocity field(LoS).
*/

 float fbox_read(float * list,const char *boxes_PATH, const char * boxname, int dim);//How to deal with array?
 void kbox_read(fftw_complex * list,const char *boxes_PATH, const char * boxname, int dim);
 double dbox_read(double * list,const char *boxes_PATH, const char * boxname, int dim);
 void fbox_write(const float * list, const char * boxname, int dim); 
 void kbox_write(const fftw_complex * list, const char * boxname, int dim); 
 void dbox_write(const double * list, const char * boxname, int dim); 
 double * dmyalloc(unsigned long long num);
 unsigned long long * ullmyalloc(unsigned long long num);
 float * fmyalloc(unsigned long long num);
 
 extern FILE *LOG;

 float fbox_read(float * list,const char *boxes_PATH, const char * boxname, int dim)//You need to allocate the memory for the the list by yourself
 {
 	FILE * fp;
 	float data_ave = 0;
 	char filename[M_BOXNAME];
 	memset(filename,'\0',M_BOXNAME);
 	sprintf(filename,"%s/%s",boxes_PATH,boxname);
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("Error: Can't open file %s\n", filename);
		fprintf(LOG, "Error: Can't open file %s\n", filename);
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	if(fread(list, sizeof(float), box_vol(dim), fp)!=(box_vol(dim))){
		printf("Error: Error happens while reading %s\n",filename);
		fprintf(LOG, "Error: Error happens while reading %s\n",filename);
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	if(fclose(fp)!=0){
		printf("Error: Can't close file %s\n",filename);
		fprintf(LOG, "Error: Can't close file %s\n",filename);
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	for(unsigned long long ct=0;ct<box_vol(dim);ct++)
		data_ave += *(list + ct);
	data_ave /= (float)box_vol(dim);
	return data_ave;
 }

 void kbox_read(fftw_complex * list,const char *boxes_PATH, const char * boxname, int dim)//You need to allocate the memory for the the list by yourself
 {
 	FILE * fp;
 	char filename[M_BOXNAME];
 	sprintf(filename,"%s%s",boxes_PATH,boxname);
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("Error: Can't open file %s\n", filename);
		fprintf(LOG, "Error: Can't open file %s\n", filename);
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	if(fread(list, sizeof(fftw_complex), kbox_vol(dim), fp)!=(kbox_vol(dim))){
		printf("Error: Error happens while reading %s\n",filename);
		fprintf(LOG, "Error: Error happens while reading %s\n",filename);
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	if(fclose(fp)!=0){
		printf("Error: Can't close file %s\n",filename);
		fprintf(LOG, "Error: Can't close file %s\n",filename);
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
 }

 double dbox_read(double * list,const char *boxes_PATH, const char * boxname, int dim)//You need to allocate the memory for the the list by yourself
 {
 	FILE * fp;
 	double data_ave = 0;
 	char filename[M_BOXNAME];
 	sprintf(filename,"%s/%s",boxes_PATH,boxname);
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("Error: Can't open file %s\n", filename);
		fprintf(LOG, "Error: Can't open file %s\n", filename);
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	if(fread(list, sizeof(double), box_vol(dim), fp)!=(box_vol(dim))){
		printf("Error: Error happens while reading %s\n",filename);
		fprintf(LOG, "Error: Error happens while reading %s\n",filename);
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	if(fclose(fp)!=0){
		printf("Error: Can't close file %s\n",filename);
		fprintf(LOG, "Error: Can't close file %s\n",filename);
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	for(unsigned long long ct=0;ct<box_vol(dim);ct++)
		data_ave += *(list + ct);
	data_ave /= (double)box_vol(dim);
	return data_ave;
 }

 void fbox_write(const float * list,const char * boxname, int dim)
 {
 	FILE * fp;
 	char filename[M_BOXNAME];
 	sprintf(filename,"%s/%s",MMRRM_BOX_OP,boxname);
 	if((fp=fopen(filename,"w"))==NULL)
 	{
 		printf("Error: Can't open file %s\n", filename);
 		fprintf(LOG, "Error: Can't open file %s\n", filename);
 		fclose(LOG);
 		exit(EXIT_FAILURE);
 	}
 	if((fwrite(list, sizeof(float), box_vol(dim), fp))!=box_vol(dim))
 	{
 		printf("Error: Can't write the output box %s\n", filename);
 		fprintf(LOG, "Error: Can't write the output box %s\n", filename);
 		fclose(LOG);
 		exit(EXIT_FAILURE);
 	}
 	if(fclose(fp)!=0){
 		printf("Error: Can't close file %s\n", filename);
 		fprintf(LOG, "Error: Can't close file %s\n", filename);
 		fclose(LOG);
 		exit(EXIT_FAILURE);
 	}
 }

 void kbox_write(const fftw_complex * list,const char * boxname, int dim)
 {
 	FILE * fp;
 	char filename[M_BOXNAME];
 	sprintf(filename,"%s/%s",MMRRM_BOX_OP,boxname);
 	if((fp=fopen(filename,"w"))==NULL)
 	{
 		printf("Error: Can't open file %s\n", filename);
 		fprintf(LOG, "Error: Can't open file %s\n", filename);
 		fclose(LOG);
 		exit(EXIT_FAILURE);
 	}
 	if((fwrite(list, sizeof(fftw_complex), kbox_vol(dim), fp))!=kbox_vol(dim))
 	{
 		printf("Error: Can't write the output box %s\n", filename);
 		fprintf(LOG, "Error: Can't write the output box %s\n", filename);
 		fclose(LOG);
 		exit(EXIT_FAILURE);
 	}
 	if(fclose(fp)!=0){
 		printf("Error: Can't close file %s\n", filename);
 		fprintf(LOG, "Error: Can't close file %s\n", filename);
 		fclose(LOG);
 		exit(EXIT_FAILURE);
 	}
 }

 void dbox_write(const double * list,const char * boxname, int dim)
 {
 	FILE * fp;
 	char filename[M_BOXNAME];
 	sprintf(filename,"%s/%s",MMRRM_BOX_OP,boxname);
 	if((fp=fopen(filename,"w"))==NULL)
 	{
 		printf("Error: Can't open file %s\n", filename);
 		fprintf(LOG, "Error: Can't open file %s\n", filename);
 		fclose(LOG);
 		exit(EXIT_FAILURE);
 	}
 	if((fwrite(list, sizeof(double), box_vol(dim), fp))!=box_vol(dim))
 	{
 		printf("Error: Can't write the output box %s\n", filename);
 		fprintf(LOG, "Error: Can't write the output box %s\n", filename);
 		fclose(LOG);
 		exit(EXIT_FAILURE);
 	}
 	if(fclose(fp)!=0){
 		printf("Error: Can't close file %s\n", filename);
 		fprintf(LOG, "Error: Can't close file %s\n", filename);
 		fclose(LOG);
 		exit(EXIT_FAILURE);
 	}
 }

double * dmyalloc(unsigned long long num){
	double *list;
	if((list = (double *)malloc(sizeof(double)*num))==NULL){
		printf("Error: Can't allocate memory\n");
		fprintf(LOG,"Error: Can't allocate memory\n");
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	for(unsigned long long ct=0;ct<num;ct++)
		list[ct] = 0.0;
	return list;
}
unsigned long long * ullmyalloc(unsigned long long num){
	unsigned long long *list;
	if((list = (unsigned long long *)malloc(sizeof(unsigned long long)*num))==NULL){
		printf("Error: Can't allocate memory\n");
		fprintf(LOG,"Error: Can't allocate memory\n");
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	for(unsigned long long ct=0;ct<num;ct++)
		list[ct] = 0;
	return list;
}

float * fmyalloc(unsigned long long num){
	float *list;
	if((list = (float *)malloc(sizeof(float)*num))==NULL){
		printf("Error: Can't allocate memory\n");
		fprintf(LOG,"Error: Can't allocate memory\n");
		fclose(LOG);
		exit(EXIT_FAILURE);
	}
	for(unsigned long long ct=0;ct<num;ct++)
		*(list + ct) = 0.0;
	return list;
}