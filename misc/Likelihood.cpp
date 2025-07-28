#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <iostream>
#include <math.h>

#include "M3Load.h"
#include "../bfgs.h"

#define NUMROWS 30
#define KCOL 21
#define NCOL 22
#define XCOL 24

void mu_phi_to_alpha_beta(double mu, double phi, double* alpha, double* beta) {
    *alpha = mu * phi;
    *beta = phi * (1 - mu);
}

void get_priors(double* priors) {
    // Generate priors
    double A_mu = 0.1, A_phi = 10;
    double q_mu = 0.2, q_phi = 5;
    double c_mu = 0.1, c_phi = 10;
    double phi_min = 2, phi_scale = 1000;

    mu_phi_to_alpha_beta(A_mu, A_phi, &priors[0], &priors[1]);
    mu_phi_to_alpha_beta(q_mu, q_phi, &priors[2], &priors[3]);
    mu_phi_to_alpha_beta(c_mu, c_phi, &priors[4], &priors[5]);
    
    priors[6] = phi_min;
    priors[7] = phi_scale;
}

double log_beta_func(double alpha, double beta){
    //fprintf(stderr,"BETAFUNC\n");
    return lgamma(alpha)+lgamma(beta)-lgamma(alpha+beta);
}

double log_beta(double x, double alpha, double beta) {

    double term1 = (beta-1)*log(-x+1) ;
    double term2 = (alpha-1)*log(x);
    double term3 = log_beta_func(alpha,beta); 

    return term1+term2-term3;
}

double log_exponential(double x, double loc, double scale){
    if (x < loc) {
        return -INFINITY;
    }
    return -(x - loc) / scale - log(scale);
}

double compute_log_likelihood(const double DMGparam[],const void *dats){
  (void) dats;
    //Compute the log-likelihood across all positions and return the result
    /*double part1=0;
    double part2=0;
    double Dx;
    double alpha;
    double beta;*/

    // DMGparam =  A q c phi
    double A = DMGparam[0]; 
    double q = DMGparam[1];
    double c = DMGparam[2];
    double phi = DMGparam[3];

    double Dx;
    double alpha;
    double beta;
    double part1;
    double part2;

    double like_sum = 0;
    for (int i = 0; i < NUMROWS; i++) {
        Dx = A * pow((1 - q), fabs(M3[i][XCOL]) - 1) + c;
        //fprintf(stderr,"DX %f \n",Dx);
        alpha = Dx * phi;
        beta = (1 - Dx) * phi;
        //double likelihood = compute_log_likelihood(A, q, c, phi, M3[i][x_col], M3[i][k_col], M3[i][N_col]);
        part1 = lgamma(M3[i][NCOL]+1)+lgamma(M3[i][KCOL]+alpha)+lgamma(M3[i][NCOL]-M3[i][KCOL]+beta)+lgamma(alpha+beta);
        part2 = lgamma(M3[i][KCOL]+1)+lgamma(M3[i][NCOL]-M3[i][KCOL]+1)+lgamma(alpha)+lgamma(beta)+lgamma(M3[i][NCOL]+alpha+beta);
        //fprintf(stderr,"XCAL %f \t M3[i][NCOL] %f \n",fabs(M3[i][KCOL]),M3[i][NCOL]);
        //fprintf(stderr,"part1 is %f \t part2 %f \n",part1,part2);
        like_sum = like_sum + (part1-part2); //(part1-part2) -> likelihood
    }
    //fprintf(stderr,"A: %0.10f \t q %0.10f \t c %0.10f \t phi %0.10f\n",A,q,c,phi);
    //fprintf(stderr,"Compute log-likelihood is %f \n",(-1)*like_sum);
    return (-1)*like_sum;
}

double compute_log_prior(double* Combined){
    //fprintf(stderr,"Calling compute_log_prior \n");
    // DMGparam =  A q c phi
    double A = Combined[0]; 
    double q = Combined[1];
    double c = Combined[2];
    double phi = Combined[3];

    // Combined 1 -> 3 is the initial estimate 
    // Combined 4 -> 11 is the priors generated from get_prior function
    double lp =log_beta(A,Combined[4],Combined[5])+log_beta(q,Combined[6],Combined[7])+log_beta(c,Combined[8],Combined[9])+log_exponential(phi,Combined[10],Combined[11]);
    //fprintf(stderr,"-------------- \n compute_log_prior %f \t %f \t %f \n --------------\n",A,priors[0],priors[1]);
    return -lp;
}

double compute_log_posterior(double* CombinedParam){
    // posterior = prior * likelihood -> log posterior = log prior + log likelihood
    //fprintf(stderr,"Calling compute_log_posterior \n");
    double log_prior = compute_log_prior(CombinedParam);
    double log_likelihood = compute_log_likelihood(CombinedParam,NULL);
    return log_likelihood+log_prior;
}

void get_stat_val(double* stats,const double Combined[]) {
    // Amu Aphi qmu qphi cmu cphi phimin phiscale
    stats[0] = Combined[0]*2;
    stats[1] = Combined[1]*3.2;
    stats[2] = Combined[2]*1.1;
    stats[3] = Combined[3]*1.3;
}

#ifdef __WITH_MAIN__

int main() {
    /*double M3[MAX_ROWS][MAX_COLS];
    char tax_id[MAX_ROWS][MAX_COLS];
    char dir[MAX_ROWS][MAX_COLS];*/

    int num_rows, num_cols;
    read_count_matrix("MycoBactBamSEOutSortMDSortN.mismatches.txt.gz", M3,tax_id,dir,&num_rows, &num_cols);
    Alter_count_matrix(M3,tax_id,dir,num_rows,num_cols);

    int n_priors=8;
    double* priors = (double*) malloc(n_priors*sizeof(double));
    get_priors(priors);
    int numpars = 4;

    double lowbound[] = {0.00000001,0.00000001,0.00000001,2};
    double upbound[] = {1-0.00000001,1-0.00000001,1-0.00000001,100000};
    int nbd[] = {2,2,2,1}; //2 is both lower/upper bound
    int noisy = 0;
    double* DMGparam = (double*) malloc(numpars*sizeof(double)); // A q c phi
    double* CombinedParam = (double*) malloc((numpars + n_priors) * sizeof(double));
    DMGparam[0] = 0.1;
    DMGparam[1] = 0.1;
    DMGparam[2] = 0.01; 
    DMGparam[3] = 1000; 
    memcpy(CombinedParam, DMGparam, numpars * sizeof(double));
    memcpy(CombinedParam + numpars, priors, n_priors * sizeof(double));

    findmax_bfgs(numpars, DMGparam, nullptr,compute_log_likelihood,nullptr,lowbound,upbound,nbd,noisy);
    fprintf(stderr," Optimized DMGparam A:%f \t q:%f \t c:%f \t phi:%f \n",DMGparam[0],DMGparam[1],DMGparam[2],DMGparam[3]);
    
    fprintf(stderr,"\n Compute LOGLIKELIHOOD OF OPTIMZED VALUES %f \n",compute_log_likelihood(DMGparam,nullptr));

    free(priors);
    free(DMGparam);

    return 0;
}
#endif

//g++ Likelihood.cpp M3Load.o bfgs.o -std=c++11 -lm -lz -D __WITH_MAIN__ -o Like  
