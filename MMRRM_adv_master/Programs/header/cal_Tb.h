// cal_Tb.h 
// head file for functions related to brightness temperature calculation
// Edit history
// 2018.4.19: Add more comments, Jiachuan Xu

#include "mmrrm_adv.h"
#include <string.h>

/********************	Function Declaration	********************/

// cal_Tb: calculate $\deltaT_b$ for a data cube
// Inputs:
//		delta:	mass overdensity, dim^3 float array, [dimensionless]
//		nHI:	neutral fraction, dim^3 float array, [dimensionless]<---------------check 
//		Ts:		spin temperature, dim^3 float array, [K]
//		vlos:	peculiar velocity along LoS, dim^3 float array, [cMpc/s]
//		LoS_state: direction of LoS, 1-x 2-y 3-z, [dimensionless]
//		zrl:	redshift, [dimensionless]
// Outputs:
//		Tb:		$\deltaT_b$, dim^3 double array, [mK or K, see definition]
//		TB_AVE:	<Tb> of the dim^3 cube, [mK or K, see definition]
// Return:
//		op_thick_c: counts of optical thick cells, float
float cal_Tb(const float *delta, const float *nHI, const float *Ts,
	const float *vlos, double *Tb, int LoS_state, float zrl, 
	double *TB_AVE);

// cal_tau: calculate optical depth (tau) for one cell 
// Inputs:
//		i_p,j_p,k_p: indices for cell, i_p-x j_p-y k_p-z
//		delta:	mass overdensity, dim^3 float array, [dimensionless]
//		nHI:	neutral fraction, dim^3 float array, [dimensionless]<---------------check 
//		Ts:		spin temperature, dim^3 float array, [K]
//		det:	the determinent |1+dv/dr/aH|(Yi Mao et.al. 2012)
//		zrl:	redshift, [dimensionless]
// Return:
//		tau:	optical depth, double, [dimensionless]
double cal_tau(int i_p, int j_p, int k_p, const float *delta, 
	const float *nHI, const float *Ts, double det, float zrl);

// swap:	switch the position of two number
// Inputs:
//		a: 		number a, double
//		b: 		number b, double
void swap(double *a, double *b);

// mo:	calculate modulo of i over dim (periodic boundary condition)
// Inputs:
//		i:		index of cell 
// Return:
//		i':		i%dim_nbody
int mo(int i);

// cal_TbV_thin: calculate $\deltaT_b*V_{cell}$ for one cell, under 
//               optical-thin approximation, used for mesh2mesh mapping
// Inputs:
//		i_p,j_p,k_p: indices for cell, i_p-x j_p-y k_p-z
//		delta:	mass overdensity, dim^3 float array, [dimensionless]
//		nHI:	neutral fraction, dim^3 float array, [dimensionless]<---------------check 
//		Ts:		spin temperature, dim^3 float array, [K]
//		zrl:	redshift, [dimensionless]
// Return:
//		TbV:	$\deltaT_b*V_{cell}$, [mK or K, see definition]
//              (NOTE: units of volumn has been normalized)
double cal_TbV_thin(int i_p, int j_p, int k_p, const float *delta, 
	const float *nHI, const float *Ts, float zrl);

// cal_TbV_thick: calculate $\deltaT_b*V_{cell}$ for one cell, without 
//               optical-thin approximation, used for mesh2mesh mapping
// Inputs:
//		i_p,j_p,k_p: indices for cell, i_p-x j_p-y k_p-z
//		Ts:		spin temperature, dim^3 float array, [K]
// 		det:	the determinent |1+dv/dr/aH|(Yi Mao et.al. 2012)
//		tau:	optical depth, [dimensionless]
//		zrl:	redshift, [dimensionless]
// Return:
//		TbV:	$\deltaT_b*V_{cell}$, [mK or K, see definition]
//              (NOTE: units of volumn has been normalized)
double cal_TbV_thick(int i_p, int j_p, int k_p, const float *Ts, 
	double det, double tau, float zrl);

///////////////////////////////////////////////////////////////////////<------------obsolete?
// cal_Tb_thin: calculate $\deltaT_b$ for one cell, under optical-thin
//              approximation
// Inputs:
//		i_p,j_p,k_p: indices for cell, i_p-x j_p-y k_p-z
//		delta:	mass overdensity, dim^3 float array, [dimensionless]
//		nHI:	neutral fraction, dim^3 float array, [dimensionless]<---------------check 
//		Ts:		spin temperature, dim^3 float array, [K]
//		zrl:	redshift, [dimensionless]
//		det:	the determinent |1+dv/dr/aH|(Yi Mao et.al. 2012)
// Return:
//		Tb:		$\deltaT_b$, [mK or K, see definition]
double cal_Tb_thin(int i_p, int j_p, int k_p, const float *delta, 
	const float *nHI, const float *Ts, float zrl, double det);

///////////////////////////////////////////////////////////////////////<------------obsolete?
// cal_Tb_thick: calculate $\deltaT_b$ for one cell, without 
//               optical-thin approximation
// Inputs:
//		i_p,j_p,k_p: indices for cell, i_p-x j_p-y k_p-z
//		Ts:		spin temperature, dim^3 float array, [K]
//		tau:	optical depth, [dimensionless]
//		zrl:	redshift, [dimensionless]
// Return:
//		Tb:		$\deltaT_b$, [mK or K, see definition]
double cal_Tb_thick(int i_p, int j_p, int k_p, const float *Ts, 
	double tau, float zrl);

// Mass_Ave: calculate average neutral fraction <xHI> weighted 
//           by cell mass
// Inputs:
//		rho:	mass overdensity, dim^3 float array, [dimensionless]
//		xHI:	neutral fraction, dim^3 float array, [dimensionless] 
//		dim_nbody: dimension of the N-body simulation data cube
//		dim_rt:	dimension of the radiation transfer data cube
// Return:
//		xHI_mass_ave: <xHI>_m, [dimensionless]
float Mass_Ave(const float * rho, const float * xHI, int dim_nbody, 
	int dim_rt);

// Eta_Ave: calculate average $\eta=(1-T_{CMB}/T_s)$
// Inputs:
//		Ts:		spin temperature, dim^3 float array, [K]
//		zrl:	redshift, [dimensionless]
//		dim: 	dimension of the array
// Return:
//		<$\eta$>: [dimensionless]
double Eta_Ave(const float * Ts, float zrl, int dim);

// Tb_hat: calculate $\hat{T_b}$, numerical rules
// Inputs:
//		zrl:	redshift, [dimensionless]
//		xHI_m_ave:	mass averaged neutral fraction, [dimensionless]
// Return:
//		$\hat{T_b}$: [mK]
double Tb_hat(float zrl, float xHI_m_ave);

/*******************	Variable Declaration	*******************/

extern int dim_nbody,dim_rt; // dimension of data cube, N-body or RT
extern float box_size; // edge length of data cube, cMpc
extern FILE *LOG; // pointer -> log file 
int dim_ratio; // dim_nbody/dim_rt, must be integer

/**********************	Function Definition	**********************/

float cal_Tb(const float * delta,const float * nHI,const float * Ts,
	const float * vlos,double * Tb,int LoS_state,float zrl, 
	double *TB_AVE)
{
	// Variable declaration & initialization
	double Tb_ave = 0;
	double len = (double)box_size/(double)dim_nbody; // cell length
	unsigned long long op_thick_c=0;	//Count for optical thick cells
	int i,j,k;
	int DIVERG_CELL = 0; // cell whose 1/determinent diverge
	dim_ratio = (int)(dim_nbody / dim_rt);
	char char_op_thick[]="Optical thick cells:";
	char char_infinite[]="Diverged cells:";

	double *tau_array; // optical depth data cube <---------------------------------Debug
	double *det_array; // determinent data cube <-----------------------------------Debug
	if(DEBUG==1){
		tau_array = (double *)malloc(sizeof(double)*box_vol(dim_nbody));
		det_array = (double *)malloc(sizeof(double)*box_vol(dim_nbody));
	}

// Start OMP block
#pragma omp parallel num_threads(nthreads) default(none) \
shared(dim_nbody,LoS_state,DIVERG_CELL,zrl,len,vlos,delta, nHI, Ts,Tb, tau_array, det_array) \
private(i,j,k) \
reduction(+:op_thick_c,Tb_ave)
{
	// omp slice prependicular to x-axis
	// I/O error could happen if LoS=x 
	// NEED ATOMIC WHILE CONDUCTING MESH2MESH MAPPING 
	# pragma omp for 
	// Start for: x 
	for(i=0;i<dim_nbody;i++){
		// Start for: y 
		for(j=0;j<dim_nbody;j++){
			// Start for: z 
			for(k=0;k<dim_nbody;k++){

				int itmp,igx,igy,igz,i_start,i_end,ig;
				double vl,vr,left,right,det,tau_c,TbV;

				// Start switch block: match LoS
				// decide LoS for each cell <---------------------------------------optimize?(No) 
				switch(LoS_state){
					case 1: // LoS=x, calculate cell x-boundary vx
					itmp=mo(i-1); // left neighbour along x 
					// left boundary velocity, [dimensionless]
					vl = ((double)vlos[in_tr(i,j,k,dim_nbody)] + 
						  (double)vlos[in_tr(itmp,j,k,dim_nbody)]) 
					     * Mpc_over_km / (2.0 * Hz(zrl) * len);
					itmp=mo(i+1); // right neighbour along x 
					// right boundary velocity, [dimensionless]
					vr = ((double)vlos[in_tr(i,j,k,dim_nbody)] + 
						  (double)vlos[in_tr(itmp,j,k,dim_nbody)]) 
					     * Mpc_over_km / (2.0 * Hz(zrl) * len);
					// mapping: s = r + v/aH
					left = (double)i+vl; // mapped left boundary
					right = (double)i+1.0+vr; // mapped right boundary 
					igy=j;
					igz=k;
					break;
					case 2: // LoS=y, calculate cell y-boundary vy 
					itmp=mo(j-1); // left neighbour along y 
					vl = ((double)vlos[in_tr(i,j,k,dim_nbody)] + 
						  (double)vlos[in_tr(i,itmp,k,dim_nbody)]) 
					     * Mpc_over_km / (2.0 * Hz(zrl) * len);
					itmp=mo(j+1); // right neighbour along y 
					vr = ((double)vlos[in_tr(i,j,k,dim_nbody)] + 
						  (double)vlos[in_tr(i,itmp,k,dim_nbody)]) 
					     * Mpc_over_km / (2.0 * Hz(zrl) * len);
					// mapping: s = r + v/aH
					left = (double)j+vl;
					right = (double)j+1.0+vr;
					igx=i;
					igz=k;
					break;
					case 3: // LoS=z, calculate cell z-boundary vz 
					itmp=mo(k-1); 
					vl = ((double)vlos[in_tr(i,j,k,dim_nbody)] + 
						  (double)vlos[in_tr(i,j,itmp,dim_nbody)]) 
					     * Mpc_over_km / (2.0 * Hz(zrl) * len);
					itmp=mo(k+1); 
					vr = ((double)vlos[in_tr(i,j,k,dim_nbody)] + 
						  (double)vlos[in_tr(i,j,itmp,dim_nbody)]) 
					     * Mpc_over_km / (2.0 * Hz(zrl) * len);
					left = (double)k+vl;
					right = (double)k+1.0+vr;
					igx=i;
					igy=j;
					break;
					default: 
					exit(EXIT_FAILURE);
				}// End switch block
				// switch boundary if needed
				if(left>right)
					swap(&left, &right);
				// the biggest integer < left boundary 
				i_start = (int)floor(left);
				// the biggest integer < right boundary
				i_end = (int)floor(right);
				// determinent |1+dv/dr/aH|
				det = right - left;
				// calculate optical depth for this cell 
				tau_c = cal_tau(i,j,k, delta, nHI, Ts, det,zrl);
				if(DEBUG==1){
					tau_array[in_tr(i,j,k,dim_nbody)] = tau_c;//<-----------------------Debug
					det_array[in_tr(i,j,k,dim_nbody)] = det;  //<-----------------------Debug
				}
				// decide whether this cell diverge 
				if(!isfinite(tau_c)){
					# pragma omp atomic
					DIVERG_CELL ++;
				}
				// If want result in redshift space,
				// conduct mesh2mesh mapping. 
				// Start redshift space block
				if(MESH2MESH){
					// Keep [1-exp(-tau)]
					if(OPTHIN==0){
						if(tau_c>CRE_TAU){
							# pragma omp atomic 
							op_thick_c++;
							TbV = cal_TbV_thick(i,j,k,Ts,det,tau_c,zrl);
						}
						else
							TbV = cal_TbV_thin(i,j,k,delta,nHI,Ts,zrl);
					}
					// Forced optical-thin
					else
						TbV = cal_TbV_thin(i,j,k,delta,nHI,Ts,zrl);
					// Start for: MESH2MESH
					for(ig=i_start;ig<=i_end;ig++){ // engine index: ig 
						// engage the clutch
						switch(LoS_state){
							case 1: // "x-gear"
							igx=mo(ig);
							break;
							case 2: // "y-gear"
							igy=mo(ig);
							break;
							case 3: // "z-gear"
							igz=mo(ig);
							break;
							default:
							exit(EXIT_FAILURE);
						}
						// for each step, decide:
						// case 1: two boundaries fall in the same cell
						if(i_start==i_end){
						#pragma omp atomic
							Tb[in_tr(igx,igy,igz,dim_nbody)] += TbV;
						}
						// case 2: at the beginning 
						else if(ig==i_start){
						#pragma omp atomic
							Tb[in_tr(igx,igy,igz,dim_nbody)] += TbV/det*((double)i_start+1.0-left);
						}
						// case 3: at the end 
						else if(ig==i_end){
						#pragma omp atomic
							Tb[in_tr(igx,igy,igz,dim_nbody)] += TbV/det*(right - (double)i_end);
						}
						// case 4: at the middle
						else{
						#pragma omp atomic
							Tb[in_tr(igx,igy,igz,dim_nbody)] += TbV/det;
						}
					}// End for: MESH2MESH
				}// End redshift space block
				// If want result in real space, 
				// skip mesh2mesh mapping, only consider modification
				// on $\tau$ and $\deltaT_b$.
				// Start real space block
				else{
					if(OPTHIN==0){
						if(tau_c>CRE_TAU){
							# pragma omp atomic
							op_thick_c++;
							TbV = cal_Tb_thick(i,j,k,Ts,tau_c,zrl);
						}
						else
							TbV = cal_Tb_thin(i, j, k, delta, nHI, Ts, zrl, det);
					}
					// Forced optical-thin
					else
						TbV = cal_Tb_thin(i, j, k, delta, nHI, Ts, zrl, det);
					Tb[in_tr(i,j,k,dim_nbody)] += TbV;
				}// End real space block
				Tb_ave += TbV;
			}// End for: z 
		}// End for: y 
	}// End for: x 
}// End OMP block

	printf("%22s %6llu\n",char_op_thick,op_thick_c);
	fprintf(LOG, "%22s %6llu\n",char_op_thick,op_thick_c);
	printf("%22s %6d\n", char_infinite, DIVERG_CELL);
	fprintf(LOG, "%22s %6d\n", char_infinite, DIVERG_CELL);
	// Output mean $\deltaT_b$
	Tb_ave /= pow((double)dim_nbody,3.0);
	*TB_AVE = Tb_ave;

	if(DEBUG==1)
	{
		// Output tau data cube <-------------------------------------------------------Debug
		FILE *fp_dbg;
		char filename_dbg[M_BOXNAME] = {'\0'};
		char flag_opthin, flag_hights, flag_space[3];
		char cmd[M_CMD], dir[M_PATH];
		if(OPTHIN)
			flag_opthin = 'T';// optical-Thin 
		else
			flag_opthin = 'F';// Full optical depth
		if(HIGHTS)
			flag_hights = 'H';// High spin temperature limit
		else
			flag_hights = 'N';// Normal spin temperature treatment
		if(MESH2MESH)
			strcpy(flag_space, "RS");// RedShift space
		else
			strcpy(flag_space, "RE");// REal space 

		sprintf(dir, "%s/tau_binary", MMRRM_BOX_OP);
		if(access(dir,F_OK)==-1){
		sprintf(cmd,"mkdir %s",dir);
		system(cmd);
		}
		sprintf(filename_dbg, 
			"%s/tau_Aprox%c%c_%sSpace_z%05.2f_dim%04d_size%04.0f_los%d",
			dir, flag_opthin,flag_hights,flag_space,
			zrl,dim_nbody,box_size,LoS_state); 
		fp_dbg = fopen(filename_dbg,"w");
		fwrite(tau_array, sizeof(double), box_vol(dim_nbody), fp_dbg);
		fclose(fp_dbg);
		free(tau_array);
		// output determinent data cube <-----------------------------------------------Debug
		sprintf(dir, "%s/det_binary", MMRRM_BOX_OP);
		if(access(dir,F_OK)==-1){
		sprintf(cmd,"mkdir %s",dir);
		system(cmd);
		}
		sprintf(filename_dbg, 
			"%s/det_Aprox%c%c_%sSpace_z%05.2f_dim%04d_size%04.0f_los%d",
			dir, flag_opthin,flag_hights,flag_space,
			zrl,dim_nbody,box_size,LoS_state); 
		fp_dbg = fopen(filename_dbg,"w");
		fwrite(det_array, sizeof(double), box_vol(dim_nbody), fp_dbg);
		fclose(fp_dbg);
		free(det_array);
	}
	// return optical-thick cell counts 
	return (float)op_thick_c;
}

double cal_tau(int i_p,int j_p,int k_p,
const float * delta,const float * nHI, const float * Ts,double det,float zrl)
{
	double tautemp;
	int i_pc = (int)(i_p / dim_ratio);// index in RT data cube
	int j_pc = (int)(j_p / dim_ratio);
	int k_pc = (int)(k_p / dim_ratio);
	double delta_d = (double)delta[in_tr(i_p,j_p,k_p,dim_nbody)];
	double nHI_d = (double)nHI[in_tr(i_pc,j_pc,k_pc,dim_rt)];	//1/cm^3<-----------???check
	double Ts_d = (double)Ts[in_tr(i_p,j_p,k_p,dim_nbody)];
	// numerical rule for optical depth
	//		    3*c^3*A10*T_star*nHI
	// tau = --------------------------, nHI comoving or physical?<-----------------???check
	//		  32*pi*nu_0^3*Ts*H(z)*det
	tautemp = 0.02388 * (Ob*h/0.02) * 
	sqrt(0.15 / Om * pow((double)(1.0 + (double)zrl),3.0) / 10.0) 
	* (1.0 + delta_d) * nHI_d / (Ts_d * det);//<------------------------------------verify this, should use constant & parameter
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

double cal_TbV_thin(int i_p,int j_p,int k_p,
	const float * delta, const float *nHI,const float * Ts,float zrl)
{
	int i_pc = (int)(i_p / dim_ratio);// index in RT data cube
	int j_pc = (int)(j_p / dim_ratio);
	int k_pc = (int)(k_p / dim_ratio);
	double TbVtemp;
	double delta_d = (double)delta[in_tr(i_p,j_p,k_p,dim_nbody)];
	double nHI_d = (double)nHI[in_tr(i_pc,j_pc,k_pc,dim_rt)];	//1/cm^3<-----------check
	double Ts_d;//maybe use eta here is better?<------------------------------------check 
	if(HIGHTS==0){
		Ts_d = (double)Ts[in_tr(i_p,j_p,k_p,dim_nbody)];
		//TbVtemp = 0.02388*(Ob*h/0.02)*sqrt(0.15/Om*(1.0+(double)zrl)/10.0)*(1.0+delta_d)*nHI_d*(1.0-T_cmb(zrl)/Ts_d);// K
		TbVtemp = 23.88*(Ob*h/0.02)*sqrt(0.15/Om*(1.0+(double)zrl)/10.0)*(1.0+delta_d)*nHI_d*(1.0-T_cmb(zrl)/Ts_d);// mK
	}
	else
		//TbVtemp = 0.02388*(Ob*h/0.02)*sqrt(0.15/Om*(1.0+(double)zrl)/10.0)*(1.0+delta_d)*nHI_d;// K
		TbVtemp = 23.88*(Ob*h/0.02)*sqrt(0.15/Om*(1.0+(double)zrl)/10.0)*(1.0+delta_d)*nHI_d;// mK
	return TbVtemp;	//mK
}


double cal_TbV_thick(int i_p,int j_p,int k_p,
	const float * Ts,double det,double tau,float zrl)
{
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
	return TbVtemp;	//mK
}

double cal_Tb_thin(int i_p,int j_p,int k_p, 
const float *delta, const float *nHI, const float *Ts, float zrl, double det)
{
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
	return Tbtemp;	//mK
}

double cal_Tb_thick(int i_p,int j_p,int k_p, 
	const float *Ts, double tau,float zrl)
{
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
	return Tbtemp;	//mK
}

float Mass_Ave(const float * rho, const float * xHI, int dim_nbody, 
	int dim_rt)
{
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

double Eta_Ave(const float * Ts, float zrl, int dim)
{
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