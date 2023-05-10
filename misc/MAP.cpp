#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <iostream>
#include <math.h>

#include "bfgs.h"
#include "M3Load.h"
#include "Likelihood.h"


void MAP(double M3[MAX_ROWS][MAX_COLS],double* MAP){
    double lowbound[] = {0.00000001,0.00000001,0.00000001,2};
    double upbound[] = {1-0.00000001,1-0.00000001,1-0.00000001,100000};
    int nbd[] = {2,2,2,1};
    int noisy = 0;
    int numpars = 4;//sizeof(MAP) / sizeof(double); // determine length of MAP array;
    double* DMGparam = (double*) malloc(numpars*sizeof(double)); // A q c phi

    DMGparam[0] = 0.1;
    DMGparam[1] = 0.1;
    DMGparam[2] = 0.01; 
    DMGparam[3] = 1000; 

    findmax_bfgs(numpars, DMGparam, nullptr,compute_log_likelihood,nullptr,lowbound,upbound,nbd,noisy);
    
    MAP[0] = DMGparam[0];
    MAP[1] = DMGparam[1];
    MAP[2] = DMGparam[2];
    MAP[3] = DMGparam[3];
    MAP[4] = compute_log_likelihood(DMGparam,nullptr);

    free(DMGparam);
}


#ifdef __WITH_MAIN__

int main() {
    /*double M3[MAX_ROWS][MAX_COLS];
    char tax_id[MAX_ROWS][MAX_COLS];
    char dir[MAX_ROWS][MAX_COLS];*/

    int num_rows, num_cols;
    read_count_matrix("MycoBactBamSEOutSortMDSortN.mismatches.txt.gz", M3,tax_id,dir,&num_rows, &num_cols);
    Alter_count_matrix(M3,tax_id,dir,num_rows,num_cols);
    
    int numpars = 5;
    double* LlhRes = (double*) malloc(numpars*sizeof(double));    
    MAP(M3,LlhRes);
    for (int i = 0; i < numpars; i++) {
        fprintf(stderr,"Value %f \n",LlhRes[i]);
    }

    return 0;
}
#endif

//g++ MAP.cpp M3Load.o Likelihood.o bfgs.o -std=c++11 -lm -lz -D __WITH_MAIN__ -o MAP 