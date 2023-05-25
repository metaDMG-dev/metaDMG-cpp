#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <iostream>
#include <math.h>
#include <cassert>

#include "bfgs.h"
#include "M3Load.h"
#include "Likelihood.h"

#define NUMROWS 30
#define KCOL 21
#define NCOL 22
#define XCOL 24
#define C_freq 19
#define G_freq 20

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

    fprintf(stderr,"I AM RUNNING findmax_bfgs NOW\n");
    findmax_bfgs(numpars, DMGparam, nullptr,compute_log_likelihood,nullptr,lowbound,upbound,nbd,noisy);
    fprintf(stderr,"OPTIM IS DONE NOW\n");

    MAP[0] = DMGparam[0];
    MAP[1] = DMGparam[1];
    MAP[2] = DMGparam[2];
    MAP[3] = DMGparam[3];
    MAP[4] = compute_log_likelihood(DMGparam,nullptr);

    free(DMGparam);
}

void MAP_stat_file(double M3[MAX_ROWS][MAX_COLS],double* LlhRes,gzFile file){

    // MAP damage, damage_std Map_damage_significance
    int N_max = M3[0][NCOL];

    // mu = A = LlhRes[0]
    //std::cout << " N MAX " << N_max << std::endl;
    double std = std::sqrt((LlhRes[0]*(1-LlhRes[0])*(LlhRes[3]+N_max))/((LlhRes[3]+1)*N_max));
    double significance = LlhRes[0]/std;

    // MAP_Damage   MAP_signifcance     MAP_damage_std
    gzprintf(file,"MAP_damage:%f \t MAP_damage_std:%f \t MAP_damage_significance:%f \t",LlhRes[0],std,significance);

    /*
    add the significance etc for the A, q, c when done
    */

    double q = LlhRes[1];
    double A = LlhRes[0];
    double c = LlhRes[2];
    double phi = LlhRes[3];
    
    // then the optimized values MAP_A MAP_q MAP_phi MAP_c 
    gzprintf(file,"A:%f \t q:%f \t phi:%f \t c:%f \t",A,q,phi,c);
    
    double alpha;
    double beta;
    double Dx_var_numerator;
    double Dx_var_deumerator;
    double Dx_std;
    double Dx_std_norm;
    for(int i = 0; i < NUMROWS;i++){
        double Dx = A * pow((1 - q), fabs(M3[i][XCOL]) - 1) + c;
        gzprintf(file,"Dx_%d:%f \t",i,Dx);
    }

    for(int i = 0; i < NUMROWS;i++){
        double Dx = A * pow((1 - q), fabs(M3[i][XCOL]) - 1) + c;
        alpha = Dx * phi;
        beta = (1 - Dx) * phi;
        Dx_var_numerator = M3[i][NCOL]*alpha*beta*(alpha+beta+M3[i][NCOL]);
        Dx_var_deumerator = pow((alpha+beta),2)*(alpha+beta+1);
        Dx_std = std::sqrt((Dx_var_numerator/Dx_var_deumerator));
        Dx_std_norm = Dx_std / M3[i][NCOL];
        gzprintf(file,"Dx_std_%d:%f \t",i,Dx_std_norm);
        //fprintf(stderr,"Dx_%d std %f \n",i,Dx_std_norm);
    }
    /*mu = A
    std = np.sqrt(A * (1 - A) * (phi + N) / ((phi + 1) * N))
    significance = mu / std
    return mu, std, significance
    
    A:  0.21290688706956776  q:  0.35714159577362814  c:  0.009605540202966379  phi:  9369.817278739834
N MAX IS  28921
 std is  0.00486586260361796
 Sig is  43.75521966264792
 
 */
}


void M3_stat_file(double M3[MAX_ROWS][MAX_COLS],gzFile file){
    
    int N_max = 0;
    
    // N_x=1_forward  N_x=1_reverse
    gzprintf(file,"Nfwd:%d \t Nrev:%d \t",M3[0][NCOL],M3[NUMROWS-1][NCOL]);

    int N_sum = 0;
    int N_sum_fwd = 0;
    int N_sum_rev = 0;
    int K_sum = 0;
    int K_sum_fwd = 0;
    int K_sum_rev = 0;
    for(int i = 0; i < NUMROWS;i++){
        N_sum += M3[i][NCOL];
        K_sum += M3[i][KCOL];
        if (i < NUMROWS/2){
            N_sum_fwd += M3[i][NCOL];
            K_sum_fwd += M3[i][KCOL];
        }
        else if (i > NUMROWS/2){
            N_sum_rev = M3[i][NCOL];
            K_sum_rev += M3[i][KCOL];
        }
    }

    gzprintf(file,"Nsum:%d \t Nfwd:%d \t Nrev:%d \t ksum:%d \t kfwd:%d \t krev:%d\t",N_sum,N_sum_fwd,N_sum_rev,K_sum,K_sum_fwd,K_sum_rev);

    // 5'
    for(int i = 0; i < NUMROWS/2;i++){
        gzprintf(file,"K_%d:%d \t",i+1,(int) M3[i][KCOL]);
    }

    for(int i = 0; i < NUMROWS/2;i++){
        gzprintf(file,"N_%d:%d \t",i+1,(int) M3[i][NCOL]);
    }

    for(int i = 0; i < NUMROWS/2;i++){
        gzprintf(file,"f_%d:%f \t",i+1, M3[i][C_freq]);
    }

    // 3'
    for(int i = NUMROWS-1; i > NUMROWS/2;i--){
        gzprintf(file,"N_%d:%d \t",i,(int) M3[i][KCOL]);
    }

     for(int i = NUMROWS-1; i > NUMROWS/2;i--){
        gzprintf(file,"N_%d:%d \t",i,(int) M3[i][NCOL]);
    }

     for(int i = NUMROWS-1; i > NUMROWS/2;i--){
        gzprintf(file,"f_%d:%f \t",i,M3[i][G_freq]);
    }
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

    fprintf(stderr,"A: %f \t q: %f \t c: %f \t phi: %f \t llh: %f \n",LlhRes[0],LlhRes[1],LlhRes[2],LlhRes[3],LlhRes[4]);

    const char* filename = "outputtest.txt.gz";
    gzFile gz = Z_NULL;
    gz = gzopen(filename,"w");
    assert(gz!=Z_NULL);
    M3_stat_file(M3,gz);
    MAP_stat_file(M3,LlhRes,gz);
    gzprintf(gz,"\n");
    gzclose(gz);
    //;
    
    /*for (int i = 0; i < numpars; i++) {
        fprintf(stderr,"Value %f \n",LlhRes[i]);
    }*/

    return 0;
}
#endif

//g++ MAP.cpp M3Load.o Likelihood.o bfgs.o -std=c++11 -lm -lz -D __WITH_MAIN__ -o MAP 