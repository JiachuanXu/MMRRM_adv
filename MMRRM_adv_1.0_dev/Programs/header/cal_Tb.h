#include "mmrrm_adv.h"
#include <string.h>
/*
	USAGE: To calculate the differential brightness temperature signal. 
*/

float cal_Tb(const float *delta,const float *nHI,const float *Ts,const float *vlos,double *Tb,int LoS_state,float zrl, double *TB_AVE);
double cal_tau(int i_p,int j_p,int k_p, const float *delta, const float *nHI, const float *Ts, double det,float zrl);
void swap(double *a, double *b);
int mo(int i); //modulo of i over dim
double cal_TbV_thin(int i_p,int j_p,int k_p, const float *delta, const float *nHI, const float *Ts,float zrl);
double cal_TbV_thick(int i_p,int j_p,int k_p, const float *Ts, double det, double tau,float zrl);
double cal_Tb_thin(int i_p,int j_p,int k_p, const float *delta, const float *nHI, const float *Ts, float zrl, double det);
double cal_Tb_thick(int i_p,int j_p,int k_p, const float *Ts, double tau,float zrl);
float Mass_Ave(const float * rho, const float * xHI, int dim_nbody, int dim_rt);
double Eta_Ave(const float * Ts, float zrl, int dim);
double Tb_hat(float zrl, float xHI_m_ave);
extern int dim_nbody,dim_rt;
extern float box_size;
int dim_ratio;
extern FILE *LOG;

float cal_Tb(const float * delta,const float * nHI,const float * Ts,const float * vlos,double * Tb,int LoS_state,float zrl, double *TB_AVE)
{
	double Tb_ave = 0;
	double len = (double)box_size/(double)dim_nbody;
	unsigned long long op_thick_c=0;	//Count for optical thick cells
	int i,j,k;
	int DIVERG_CELL = 0;
	dim_ratio = (int)(dim_nbody / dim_rt);
	char char_op_thick[]="Optical thick cells:";
	char char_infinite[]="Diverged cells:";

	// For Debug
	double *tau_array;
	tau_array = (double *)malloc(sizeof(double)*box_vol(dim_nbody));

	/********************** Looping *****************/ 

#pragma omp parallel num_threads(nthreads) default(none) shared(dim_nbody,LoS_state,DIVERG_CELL,zrl,len,vlos,delta, nHI, Ts,Tb, tau_array) private(i,j,k) reduction(+:op_thick_c,Tb_ave)
{
	# pragma omp for 
	for(i=0;i<dim_nbody;i++){
		for(j=0;j<dim_nbody;j++){
			for(k=0;k<dim_nbody;k++){
				int itmp,igx,igy,igz,i_start,i_end,ig;
				double vl,vr,left,right,det,tau_c,TbV;
				switch(LoS_state){
					case 1:
					itmp=mo(i-1); vl = ((double)vlos[in_tr(i,j,k,dim_nbody)] + (double)vlos[in_tr(itmp,j,k,dim_nbody)]) * Mpc_over_km / (2.0 * Hz(zrl) * len);
					itmp=mo(i+1); vr = ((double)vlos[in_tr(i,j,k,dim_nbody)] + (double)vlos[in_tr(itmp,j,k,dim_nbody)]) * Mpc_over_km / (2.0 * Hz(zrl) * len);
					left = (double)i+vl;
					right = (double)i+1.0+vr;
					igy=j;
					igz=k;
					break;
					case 2:
					itmp=mo(j-1); vl = ((double)vlos[in_tr(i,j,k,dim_nbody)] + (double)vlos[in_tr(i,itmp,k,dim_nbody)]) * Mpc_over_km / (2.0 * Hz(zrl) * len);
					itmp=mo(j+1); vr = ((double)vlos[in_tr(i,j,k,dim_nbody)] + (double)vlos[in_tr(i,itmp,k,dim_nbody)]) * Mpc_over_km / (2.0 * Hz(zrl) * len);
					left = (double)j+vl;
					right = (double)j+1.0+vr;
					igx=i;
					igz=k;
					break;
					case 3:
					itmp=mo(k-1); vl = ((double)vlos[in_tr(i,j,k,dim_nbody)] + (double)vlos[in_tr(i,j,itmp,dim_nbody)]) * Mpc_over_km / (2.0 * Hz(zrl) * len);
					itmp=mo(k+1); vr = ((double)vlos[in_tr(i,j,k,dim_nbody)] + (double)vlos[in_tr(i,j,itmp,dim_nbody)]) * Mpc_over_km / (2.0 * Hz(zrl) * len);
					left = (double)k+vl;
					right = (double)k+1.0+vr;
					igx=i;
					igy=j;
					break;
					default:
					exit(EXIT_FAILURE);
				}
				if(left>right)
					swap(&left, &right);
				i_start = (int)floor(left);
				i_end = (int)floor(right);
				det = right - left;

				tau_c = cal_tau(i,j,k, delta, nHI, Ts, det,zrl);
				// For Debug
				tau_array[in_tr(i,j,k,dim_nbody)] = tau_c;

				if(!isfinite(tau_c)){
					# pragma omp atomic
					DIVERG_CELL ++;
				}

				if(MESH2MESH){	// Redshift space
					if(OPTHIN==0){
						if(tau_c>CRE_TAU){
							op_thick_c++;
							TbV = cal_TbV_thick(i,j,k,Ts,det,tau_c,zrl);
						}
						else
							TbV = cal_TbV_thin(i,j,k,delta,nHI,Ts,zrl);
					}
					//Artificially adopt optical thin approximation
					else
						TbV = cal_TbV_thin(i,j,k,delta,nHI,Ts,zrl);
					/****************** Mesh to Mesh: In *******************/
					for(ig=i_start;ig<=i_end;ig++){
						switch(LoS_state){
							case 1:
							igx=mo(ig);
							break;
							case 2:
							igy=mo(ig);
							break;
							case 3:
							igz=mo(ig);
							break;
							default:
							exit(EXIT_FAILURE);
						}
						if(i_start==i_end){
						#pragma omp atomic
							Tb[in_tr(igx,igy,igz,dim_nbody)] += TbV;
						}
						else if(ig==i_start){
						#pragma omp atomic
							Tb[in_tr(igx,igy,igz,dim_nbody)] += TbV/det*((double)i_start+1.0-left);
						}
						else if(ig==i_end){
						#pragma omp atomic
							Tb[in_tr(igx,igy,igz,dim_nbody)] += TbV/det*(right - (double)i_end);
						}
						else{
						#pragma omp atomic
							Tb[in_tr(igx,igy,igz,dim_nbody)] += TbV/det;
						}
					}
				/**************** Mesh to Mesh: Out ***************/
				}
				else{	// Real Space
					if(OPTHIN==0){
						if(tau_c>CRE_TAU){
							op_thick_c++;
							TbV = cal_Tb_thick(i,j,k,Ts,tau_c,zrl);
						}
						else
							TbV = cal_Tb_thin(i, j, k, delta, nHI, Ts, zrl, det);
					}
					//Artificially adopt optical thin approximation
					else
						TbV = cal_Tb_thin(i, j, k, delta, nHI, Ts, zrl, det);
					Tb[in_tr(i,j,k,dim_nbody)] += TbV;
				}
				Tb_ave += TbV;
			}
		}
	}
}
	/************************** Looping out ***************************/

	printf("%22s %6llu\n",char_op_thick,op_thick_c);
	fprintf(LOG, "%22s %6llu\n",char_op_thick,op_thick_c);
	printf("%22s %6d\n", char_infinite, DIVERG_CELL);
	fprintf(LOG, "%22s %6d\n", char_infinite, DIVERG_CELL);
	Tb_ave /= pow((double)dim_nbody,3.0);
	*TB_AVE = Tb_ave;

	// For Debug
	FILE *fp_dbg;
	char filename_dbg[M_BOXNAME] = {'\0'};
	char flag_opthin, flag_hights, flag_space[3];
	if(OPTHIN)
		flag_opthin = 'T';
	else
		flag_opthin = 'F';
	if(HIGHTS)
		flag_hights = 'H';
	else
		flag_hights = 'N';
	if(MESH2MESH)
		strcpy(flag_space, "RS");
	else
		strcpy(flag_space, "RE");
	sprintf(filename_dbg, "../Output_boxes/tau_binary/tau_pdf_debug_Aprox%c%c_%sSpace_z%05.2f_dim%04d_size%04.0f_los%d",flag_opthin,flag_hights,flag_space,zrl,dim_nbody,box_size,LoS_state); 
	fp_dbg = fopen(filename_dbg,"w");
	fwrite(tau_array, sizeof(double), box_vol(dim_nbody), fp_dbg);
	fclose(fp_dbg);
	free(tau_array);

	return (float)op_thick_c;
}

double cal_tau(int i_p,int j_p,int k_p,const float * delta,const float * nHI,const float * Ts,double det,float zrl){
	double tautemp;
	int i_pc = (int)(i_p / dim_ratio);
	int j_pc = (int)(j_p / dim_ratio);
	int k_pc = (int)(k_p / dim_ratio);
	double delta_d = (double)delta[in_tr(i_p,j_p,k_p,dim_nbody)];
	double nHI_d = (double)nHI[in_tr(i_pc,j_pc,k_pc,dim_rt)];	//1/cm^3
	double Ts_d;
	//if(HIGHTS==0)
		Ts_d = (double)Ts[in_tr(i_p,j_p,k_p,dim_nbody)];
	//else
	//	Ts_d = Ts_lim;
	tautemp = 0.02388 * (Ob*h/0.02) * sqrt(0.15 / Om * pow((double)(1.0 + (double)zrl),3.0) / 10.0) * (1.0 + delta_d) * nHI_d / (Ts_d * det);
	return tautemp;
}

void swap(double *a, double *b){
	double t;
	t=*a;
	*a=*b;
	*b=t;
}

int mo(int i){
	extern int dim_nbody;
	while(i<0){
		i=i+dim_nbody;
	}
	i=i%dim_nbody;
	return i;
}

double cal_TbV_thin(int i_p,int j_p,int k_p,const float * delta,const float *nHI,const float * Ts,float zrl){
	int i_pc = (int)(i_p / dim_ratio);
	int j_pc = (int)(j_p / dim_ratio);
	int k_pc = (int)(k_p / dim_ratio);
	double TbVtemp;
	double delta_d = (double)delta[in_tr(i_p,j_p,k_p,dim_nbody)];
	double nHI_d = (double)nHI[in_tr(i_pc,j_pc,k_pc,dim_rt)];	//1/cm^3
	double Ts_d;
	if(HIGHTS==0){
		Ts_d = (double)Ts[in_tr(i_p,j_p,k_p,dim_nbody)];
		//TbVtemp = 0.02388*(Ob*h/0.02)*sqrt(0.15/Om*(1.0+(double)zrl)/10.0)*(1.0+delta_d)*nHI_d*(1.0-T_cmb(zrl)/Ts_d);// K
		TbVtemp = 23.88*(Ob*h/0.02)*sqrt(0.15/Om*(1.0+(double)zrl)/10.0)*(1.0+delta_d)*nHI_d*(1.0-T_cmb(zrl)/Ts_d);// mK
	}
	else
		//TbVtemp = 0.02388*(Ob*h/0.02)*sqrt(0.15/Om*(1.0+(double)zrl)/10.0)*(1.0+delta_d)*nHI_d;// K
		TbVtemp = 23.88*(Ob*h/0.02)*sqrt(0.15/Om*(1.0+(double)zrl)/10.0)*(1.0+delta_d)*nHI_d;// mK
	return TbVtemp;	//K
}


double cal_TbV_thick(int i_p,int j_p,int k_p,const float * Ts,double det,double tau,float zrl){
	double TbVtemp;
	double Ts_d;
	double mtau = 0.0 - tau;
	Ts_d = (double)Ts[in_tr(i_p,j_p,k_p,dim_nbody)];
	if(HIGHTS==0){
		TbVtemp = (1.0-exp(mtau))*(Ts_d-T_cmb(zrl))*det*1000.0/(1.0+(double)zrl);//mK
	}
	else
		TbVtemp = (1.0-exp(mtau))*(Ts_d)*det*1000.0/(1.0+(double)zrl);//mK
	//TbVtemp = (1.0-exp(mtau))*(Ts_d-T_cmb(zrl))*det/(1.0+(double)zrl);// K
	//TbVtemp = (1.0-exp(mtau))*(Ts_d-T_cmb(zrl))*det*1000.0/(1.0+(double)zrl);//mK
	return TbVtemp;	//K
}

double cal_Tb_thin(int i_p,int j_p,int k_p, const float *delta, const float *nHI, const float *Ts, float zrl, double det){
	int i_pc = (int)(i_p / dim_ratio);
	int j_pc = (int)(j_p / dim_ratio);
	int k_pc = (int)(k_p / dim_ratio);
	double Tbtemp;
	double delta_d = (double)delta[in_tr(i_p,j_p,k_p,dim_nbody)];
	double nHI_d = (double)nHI[in_tr(i_pc,j_pc,k_pc,dim_rt)];	//1/cm^3
	double Ts_d;
	if(HIGHTS==0){
		Ts_d = (double)Ts[in_tr(i_p,j_p,k_p,dim_nbody)];
		//Tbtemp = 0.02388*(Ob*h/0.02)*sqrt(0.15/Om*(1.0+(double)zrl)/10.0)*(1.0+delta_d)*nHI_d*(1.0-T_cmb(zrl)/Ts_d);// K
		Tbtemp = 23.88 * (Ob * h / 0.02) * sqrt(0.15 / Om * (1.0 + (double)zrl) / 10.0) * (1.0 + delta_d) * nHI_d * (1.0 - T_cmb(zrl) / Ts_d) / det;// mK
	}
	else
		//Tbtemp = 0.02388*(Ob*h/0.02)*sqrt(0.15/Om*(1.0+(double)zrl)/10.0)*(1.0+delta_d)*nHI_d;// K
		Tbtemp = 23.88 * (Ob * h / 0.02) * sqrt(0.15 / Om * (1.0 + (double)zrl) / 10.0) * (1.0 + delta_d) * nHI_d / det;// mK
	return Tbtemp;	//K
}

double cal_Tb_thick(int i_p,int j_p,int k_p, const float *Ts, double tau,float zrl){
	double Tbtemp;
	double Ts_d;
	double mtau = 0.0 - tau;
	Ts_d = (double)Ts[in_tr(i_p,j_p,k_p,dim_nbody)];
	if(HIGHTS==0){
		Tbtemp = (1.0 - exp(mtau)) * (Ts_d - T_cmb(zrl)) * 1000.0 / (1.0 + (double)zrl);//mK
	}
	else
		Tbtemp = (1.0 - exp(mtau)) * (Ts_d) * 1000.0 / (1.0 + (double)zrl);//mK
	//Tbtemp = (1.0-exp(mtau))*(Ts_d-T_cmb(zrl))*det/(1.0+(double)zrl);// K
	//Tbtemp = (1.0-exp(mtau))*(Ts_d-T_cmb(zrl))*det*1000.0/(1.0+(double)zrl);//mK
	return Tbtemp;	//K
}

float Mass_Ave(const float * rho, const float * xHI, int dim_nbody, int dim_rt){
	float H_Mass = 0.0;
	float HI_Mass = 0.0;
	float xHI_mass_ave;
	int dim_ratio = (int)(dim_nbody / dim_rt);
	for(unsigned long long ct=0;ct<box_vol(dim_nbody);ct++){
		H_Mass += (1.0 + *(rho + ct));
		HI_Mass += (*(rho + ct) + 1.0)*(*(xHI + (int)(ct/dim_ratio)));
	}
	xHI_mass_ave = HI_Mass / H_Mass;
	return xHI_mass_ave;
}

double Eta_Ave(const float * Ts, float zrl, int dim){
	double ave=0;
	for(unsigned long long ct = 0; ct < box_vol(dim); ct++)
		ave += (1.0 - (double)T_cmb(zrl)/(double)Ts[ct]);
	ave /= (double)box_vol(dim);
	return ave;
}

double Tb_hat(float zrl, float xHI_m_ave){
	double Th;
	Th = 0.02388*(Ob*h*h/0.02)*sqrt(0.15/(Om*h*h)*(1.0+(double)zrl)/10.0)*(double)xHI_m_ave;
	return Th;
}