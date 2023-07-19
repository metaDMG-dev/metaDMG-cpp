#include <cmath>
#include <cassert>
#include <cstring>
#include <vector>
#include <zlib.h>

#include "bfgs.h"


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


double **read1_ugly_matrix(const char *fname){
  fprintf(stderr,"\t-> Reading file: \'%s\'\n",fname);
  gzFile gz = Z_NULL;
  assert((gz=gzopen(fname,"rb"))!=Z_NULL);
  char buf[4096];
  char *taxid = NULL;
  gzgets(gz,buf,4096);
  std::vector<double> xcol;
  std::vector<double> kcol;
  std::vector<double> ncol;
  
  while(gzgets(gz,buf,4096)){
    char *tok = strtok(buf,"\t\n ");
    if(taxid==NULL)
      taxid = strdup(tok);
    else
      assert(strcmp(taxid,tok)==0);

    tok = strtok(NULL,"\t\n ");
    int is5 = 1;
    if(strcmp(tok,"3'")==0)
      is5 = 0;

    int pos = atoi(strtok(NULL,"\t\n "));
    double data[16];
    for(int i=0;i<16;i++)
      data[i] = atof(strtok(NULL,"\t\n "));

    double nsum = 0;
    double k =0;
    for(int i=0;i<4;i++)
      if(is5)
	nsum += data[i+4];
      else
	nsum += data[i+8];
    if(is5)
      k = data[7];
    else
      k = data[8];
    xcol.push_back(pos);
    ncol.push_back(nsum);
    kcol.push_back(k);
  }
  fprintf(stderr,"\t-> Number of datapoints: %lu \n",xcol.size());
  double **ret = new double*[3];
  ret[0] = new double[xcol.size()];
  ret[0][0] = xcol.size();
  ret[1] = new double[xcol.size()];
  ret[2] = new double[xcol.size()];
  for(int i=0;i<xcol.size();i++){
    ret[0][i+1] = xcol[i];
    ret[1][i] = kcol[i];
    ret[2][i] = ncol[i];
  }
  
  return ret;
}

double compute_log_likelihood(const double DMGparam[], const void *dats){

  const double **tmp =(const double **) dats;
  const double *XCOL = tmp[0];
  const double *KCOL = tmp[1];
  const double *NCOL = tmp[2];
  int NUMROWS = XCOL[0];

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
        Dx = A * pow((1 - q), fabs(XCOL[i+1]) - 1) + c;
        fprintf(stderr,"DX %f \n",Dx);
        alpha = Dx * phi;
        beta = (1 - Dx) * phi;
        //double likelihood = compute_log_likelihood(A, q, c, phi, M3[i][x_col], M3[i][k_col], M3[i][N_col]);
        part1 = lgamma(NCOL[i]+1)+lgamma(KCOL[i]+alpha)+lgamma(NCOL[i]-KCOL[i]+beta)+lgamma(alpha+beta);
        part2 = lgamma(KCOL[i]+1)+lgamma(NCOL[i]-KCOL[i]+1)+lgamma(alpha)+lgamma(beta)+lgamma(NCOL[i]+alpha+beta);
        //fprintf(stderr,"XCAL %f \t M3[i][NCOL] %f \n",fabs(M3[i][KCOL]),M3[i][NCOL]);
        //fprintf(stderr,"part1 is %f \t part2 %f \n",part1,part2);
        like_sum = like_sum + (part1-part2); //(part1-part2) -> likelihood
    }
    exit(0);
    //fprintf(stderr,"A: %0.10f \t q %0.10f \t c %0.10f \t phi %0.10f\n",A,q,c,phi);
    //fprintf(stderr,"Compute log-likelihood is %f \n",(-1)*like_sum);
    return (-1)*like_sum;
}


double optimoptim(double **dat){
 
    double priors[8];
    get_priors(priors);

    int numpars = 4;

    double lowbound[] = {0.00000001,0.00000001,0.00000001,2};
    double upbound[] = {1-0.00000001,1-0.00000001,1-0.00000001,100000};
    int nbd[] = {2,2,2,1}; //2 is both lower/upper bound
    int noisy = -1;
    
    double* DMGparam = (double*) malloc(numpars*sizeof(double)); // A q c phi
    double* CombinedParam = (double*) malloc((numpars + 8) * sizeof(double));
    DMGparam[0] = 0.1;
    DMGparam[1] = 0.1;
    DMGparam[2] = 0.01; 
    DMGparam[3] = 1000; 
    memcpy(CombinedParam, DMGparam, numpars * sizeof(double));
    memcpy(CombinedParam + numpars, priors, 8 * sizeof(double));

    double lik = findmax_bfgs(numpars, DMGparam, dat,compute_log_likelihood,nullptr,lowbound,upbound,nbd,noisy);
    fprintf(stderr,"\t-> Optimized DMGparam A:%f \t q:%f \t c:%f \t phi:%f lik:%f \n",DMGparam[0],DMGparam[1],DMGparam[2],DMGparam[3],lik);
    
    free(DMGparam);
    return lik;
}




#ifdef __WITH_MAIN__

int main(int argc,char **argv){
  const char *fname = "MycoBactBamSEOutSortMDSortN.mismatches.txt.gz";
  if(argc>1)
    fname = strdup(argv[1]);
  double **dat = read1_ugly_matrix(fname);
  optimoptim(dat);
}

#endif
