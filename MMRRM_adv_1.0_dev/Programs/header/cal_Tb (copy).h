#include "mmrrm_adv.h"

/*
	USAGE: To calculate the differential brightness temperature signal. 
*/

double cal_Tb(const float *delta,const float *nHI,const float *Ts,const float *vlos,double *Tb,int LoS_state,float zrl);
double cal_tau(int i_p,int j_p,int k_p, const float *delta, const float *nHI, const float *Ts, double det,float zrl);
void swap(double *a, double *b);
int mo(int i); //modulo of i over dim
double cal_TbV_thin(int i_p,int j_p,int k_p, const float *delta, const float *nHI, const float *Ts,float zrl);
double cal_TbV_thick(int i_p,int j_p,int k_p, const float *Ts, double det, double tau,float zrl);
extern int dim;
extern float box_size;
extern FILE *LOG;

double cal_Tb(const float * delta,const float * nHI,const float * Ts,const float * vlos,double * Tb,int LoS_state,float zrl)
{
	
	double Tb_ave=0;
	double len = (double)box_size/(double)dim;
	unsigned long long op_thick_c=0;	//Count for optical thick cells
	int i_dyn,j_dyn,k_dyn;
	int *i_shell,*j_shell,*k_shell;
	char char_op_thick[]="Optical thick cells:";
switch(LoS_state){
						case 1:// LoS=x
						i_shell = &k_dyn;
						j_shell = &i_dyn;
						k_shell = &j_dyn;
						break;
						case 2:// LoS=y
						j_shell = &k_dyn;
						k_shell = &i_dyn;
						i_shell = &j_dyn;
						break;
						case 3:// LoS=z
						k_shell = &k_dyn;
						i_shell = &i_dyn;
						j_shell = &j_dyn;
						break;
						default:
						printf("cal_Tb.h: Error: Illegal LoS_state value!");
						exit(EXIT_FAILURE);
					}

//# pragma omp parallel for num_threads(nthreads) reduction(+:op_thick_c) reduction(+:Tb_ave) private(i_dyn,j_dyn,k_dyn,i_shell,j_shell,k_shell,JOINT) shared(dim,LoS_state,zrl,len,Ts,nHI,delta,vlos,Tb)
	for(i_dyn=0;i_dyn<dim;i_dyn++){
		for(j_dyn=0;j_dyn<dim;j_dyn++){
			for(k_dyn=0;k_dyn<dim;k_dyn++){
				int i_start,i_end,itmp;
				double vl,vr,left,right,det,tau_c,TbV;
				itmp = k_dyn;
				
					
				

				vl = ((double)vlos[in_tr(mo(*i_shell),mo(*j_shell),mo(*k_shell),dim)])*Mpc_over_km/(2.0*Hz(zrl)*len);
				k_dyn--;
				vl += ((double)vlos[in_tr(mo(*i_shell),mo(*j_shell),mo(*k_shell),dim)])*Mpc_over_km/(2.0*Hz(zrl)*len);
				k_dyn++;
				vr = ((double)vlos[in_tr(mo(*i_shell),mo(*j_shell),mo(*k_shell),dim)])*Mpc_over_km/(2.0*Hz(zrl)*len);
				k_dyn++;
				vr += ((double)vlos[in_tr(mo(*i_shell),mo(*j_shell),mo(*k_shell),dim)])*Mpc_over_km/(2.0*Hz(zrl)*len);
				k_dyn--;
				left = (double)k_dyn + vl;
				right = (double)k_dyn + 1.0 + vr;
				if(left>right)
					swap(&left, &right);
				i_start = (int)floor(left);
				i_end = (int)floor(right);
				det = right - left;
				tau_c = cal_tau(*i_shell,*j_shell,*k_shell, delta, nHI, Ts, det,zrl);
				if(OPTHIN==0){
					if(tau_c > CRE_TAU){
						op_thick_c += 1;
						TbV = cal_TbV_thick(*i_shell,*j_shell,*k_shell,Ts,det,tau_c,zrl);
					}
					else
						TbV = cal_TbV_thin(*i_shell,*j_shell,*k_shell,delta,nHI,Ts,zrl);
				}
				//Artificially adopt optical thin approximation
				else
					TbV = cal_TbV_thin(*i_shell,*j_shell,*k_shell,delta,nHI,Ts,zrl);
				for(k_dyn = i_start;k_dyn <= i_end; k_dyn++){
					if(i_start==i_end)
						Tb[in_tr(mo(*i_shell),mo(*j_shell),mo(*k_shell),dim)] += TbV;
					else if(k_dyn==i_start)
						Tb[in_tr(mo(*i_shell),mo(*j_shell),mo(*k_shell),dim)] += TbV/det*((double)i_start+1.0-left);
					else if(k_dyn==i_end)
						Tb[in_tr(mo(*i_shell),mo(*j_shell),mo(*k_shell),dim)] += TbV/det*(right - (double)i_end);
					else
						Tb[in_tr(mo(*i_shell),mo(*j_shell),mo(*k_shell),dim)] += TbV/det;
				}
				k_dyn = itmp;
				Tb_ave += TbV;
			}
		}

	}
	printf("Brightness temperature generated successfully!\n");
	fprintf(LOG, "Brightness temperature generated successfully!\n");
	printf("%22s %6llu\n",char_op_thick,op_thick_c);
	fprintf(LOG, "%22s %6llu\n",char_op_thick,op_thick_c);
	Tb_ave /= pow((double)dim,3.0);
	return Tb_ave;
}

double cal_tau(int i_p,int j_p,int k_p,const float * delta,const float * nHI,const float * Ts,double det,float zrl){
	double tautemp;
	double delta_d = (double)delta[in_tr(i_p,j_p,k_p,dim)];
	double nHI_d = (double)nHI[in_tr(i_p,j_p,k_p,dim)];	//1/cm^3
	double Ts_d;
	//if(HIGHTS==0)
		Ts_d = (double)Ts[in_tr(i_p,j_p,k_p,dim)];
	//else
	//	Ts_d = Ts_lim;
	tautemp = 0.02388*(Ob*h/0.02)*sqrt(0.15/Om*pow((double)(1.0+(double)zrl),3.0)/10.0)*(1.0+delta_d)*nHI_d/(Ts_d*det);
	return tautemp;
}

void swap(double *a, double *b){
	double t;
	t=*a;
	*a=*b;
	*b=t;
}

int mo(int i){
	extern int dim;
	while(i<0){
		i=i+dim;
	}
	i=i%dim;
	return i;
}

double cal_TbV_thin(int i_p,int j_p,int k_p,const float * delta,const float *nHI,const float * Ts,float zrl){
	double TbVtemp;
	double delta_d = (double)delta[in_tr(i_p,j_p,k_p,dim)];
	double nHI_d = (double)nHI[in_tr(i_p,j_p,k_p,dim)];	//1/cm^3
	double Ts_d;
	if(HIGHTS==0){
		Ts_d = (double)Ts[in_tr(i_p,j_p,k_p,dim)];
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
	double mtau = 0.0-tau;
	Ts_d = (double)Ts[in_tr(i_p,j_p,k_p,dim)];
	if(HIGHTS==0){
		TbVtemp = (1.0-exp(mtau))*(Ts_d-T_cmb(zrl))*det*1000.0/(1.0+(double)zrl);//mK
	}
	else
		TbVtemp = (1.0-exp(mtau))*(Ts_d)*det*1000.0/(1.0+(double)zrl);//mK
	//TbVtemp = (1.0-exp(mtau))*(Ts_d-T_cmb(zrl))*det/(1.0+(double)zrl);// K
	//TbVtemp = (1.0-exp(mtau))*(Ts_d-T_cmb(zrl))*det*1000.0/(1.0+(double)zrl);//mK
	return TbVtemp;	//K
}