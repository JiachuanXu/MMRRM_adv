#include <stdio.h>
#include <stdlib.h>

int main(){
	char filename[]= "args_mmrrm_adv.dat";
	FILE *fp;
	float z= (6*1.0001+1)/1.02-1;
	fp = fopen(filename, "w");
	fprintf(fp,"~/Desktop/21cmFAST-master/Boxes\n256\n256\n300\n83");
	for(int i=0;i<83;i++){
		z=(1+z)*1.02-1;
		fprintf(fp,"\nupdated_smoothed_deltax_z%06.2f_256_300Mpc",z);
		fprintf(fp,"\nxH_nohalos_z%06.2f_256_300Mpc",z);
		fprintf(fp,"\nTs_z%06.2f_zetaX2.0e+56_alphaX1.2_TvirminX3.0e+04_zetaIon40.00_Pop2_256_300Mpc",z);
		fprintf(fp,"\nupdated_vx_z%06.2f_256_300Mpc",z);
		fprintf(fp,"\nupdated_vy_z%06.2f_256_300Mpc",z);
		fprintf(fp,"\nupdated_vz_z%06.2f_256_300Mpc",z);
		fprintf(fp,"\n%.2f",z);
	}
	fclose(fp);
	return 0;
}
