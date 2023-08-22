#include <cmath>
#include <cassert>
#include <cstring>
#include <vector>
#include <cstdio>
#include <zlib.h>
#include <iostream>
#include "dfit_optim.h"
#include "bfgs.h"

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
  ret[0] = new double[xcol.size()+2];
  ret[0][0] = xcol.size();
  ret[0][1] = 0;//<- function counter
  ret[1] = new double[xcol.size()];
  ret[2] = new double[xcol.size()];
  for(int i=0;i<xcol.size();i++){
    ret[0][i+2] = xcol[i];
    ret[1][i] = kcol[i];
    ret[2][i] = ncol[i];
    //   fprintf(stderr,"-> %d) %f %f %f\n",i,ret[0][i+1],ret[1][i],ret[2][i]);
  }
  
  return ret;
}

double compute_log_likelihood(const double DMGparam[], const void *dats){
  
  double **tmp =(double **) dats;
  const double *XCOL = tmp[0] +2;
  const double *KCOL = tmp[1];
  const double *NCOL = tmp[2];
  int NUMROWS = tmp[0][0];
  tmp[0][1] = tmp[0][1]+1;//increment llh function counter for nice information
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
        Dx = A * pow((1 - q), fabs(XCOL[i])) + c;
	      //fprintf(stderr,"DX %f \n",Dx);
        alpha = Dx * phi;
        beta = (1 - Dx) * phi;
	      //fprintf(stderr,"[%d] XCOL: %f KCOL: %f NCOL: %f\n",i,XCOL[i],KCOL[i],NCOL[i]);
        //double likelihood = compute_log_likelihood(A, q, c, phi, M3[i][x_col], M3[i][k_col], M3[i][N_col]);
        part1 = lgamma(NCOL[i]+1)+lgamma(KCOL[i]+alpha)+lgamma(NCOL[i]-KCOL[i]+beta)+lgamma(alpha+beta);
        part2 = lgamma(KCOL[i]+1)+lgamma(NCOL[i]-KCOL[i]+1)+lgamma(alpha)+lgamma(beta)+lgamma(NCOL[i]+alpha+beta);
        //fprintf(stderr,"XCAL %f \t M3[i][NCOL] %f \n",fabs(M3[i][KCOL]),M3[i][NCOL]);
	      //fprintf(stderr,"part1 is %f \t part2 %f \n",part1,part2);
        like_sum = like_sum + (part1-part2); //(part1-part2) -> likelihood
    }
    //    exit(0);
    //  fprintf(stderr,"A: %0.10f \t q %0.10f \t c %0.10f \t phi %0.10f\n",A,q,c,phi);
    //fprintf(stderr,"Compute log-likelihood is %f \n",(-1)*like_sum);
    
    return -like_sum;
}

//invec is 4 double long, pre values are start, post value are the optimized parameters
// A q c phi
double optim1(double *invec,double **dat){
  int numpars = 4;
  double lowbound[] = {0.00000001,0.00000001,0.00000001,2};
  double upbound[] = {1-0.00000001,1-0.00000001,0.25,100000};
  int nbd[] = {2,2,2,1}; //2 is both lower/upper bound
  int noisy = -1;
  //  dat[0][1] = 1;
  double lik = findmax_bfgs(numpars, invec, dat,compute_log_likelihood,nullptr,lowbound,upbound,nbd,noisy);
  // fprintf(stderr,"\t-> Optimized DMGparam A:%f \t q:%f \t c:%f \t phi:%f lik:%f nit: %.0f\n",invec[0],invec[1],invec[2],invec[3],lik,dat[0][1]);
  invec[4] = lik;
  invec[5] = dat[0][1];//<- plugin ncall of the objective function
  return lik;
}

double optimoptim(double *invec,double **dat,int nopt){
  double llhs=optim1(invec,dat);
  double pars[6];
  for(int i=0;i<nopt;i++){
    pars[0] = drand48()*0.8+0.1;
    pars[1] = drand48();
    pars[2] = drand48()*0.1+0.01;
    pars[3] = drand48()*10;
    
    double lik = optim1(pars,dat);
    if(lik>llhs){
      //      fprintf(stderr,"swapping\n");
      llhs=lik;
      memcpy(invec,pars,sizeof(double)*6);
    }
  }
    
  return llhs;
}


/*
  function to fill in
  std::significance::dxfit_{0...nrows}::dxfix_conf_{0...nrows}
  statpars is therefore of length 2+2*numrows;
 */
void getstat(double **dat,double *pars,double *statpars){
    // MAP damage, damage_std Map_damage_significance
    double A = pars[0];
    double q = pars[1];
    double c = pars[2];
    double phi = pars[3];
    double llh = pars[4];
    int NUMROWS = dat[0][0];
    double *XCOL = dat[0]+2;
    double *KCOL = dat[1];
    double *NCOL = dat[2];

    double N_pos = dat[2][0];
    if (N_pos < 1){
      N_pos = 1;
    }
    //  id	    A        	q	        c	          phi     	llh     	ncall	sigmaD	Zfit
    // 144905	1.000000	0.000000	0.250000	1000.000000	0.223144	27.000000	inf	0.000000
    /*fprintf(stderr,"A %f \t %f \t %f\n",A,phi,N_pos);
    fprintf(stderr,"Num 1 %f \n",A*(1-A));
    fprintf(stderr,"Num 2 %f \n",phi+N_pos);
    fprintf(stderr,"Num 2 %f \n",dat[2][0]);
    fprintf(stderr,"Num 3 %f \n",A*(1-A)*(phi+N_pos));
    fprintf(stderr,"Num 4 %f \n",(phi+1));
    fprintf(stderr,"Num 4 %f \n",(phi+1)*N_pos);*/
    double std = std::sqrt((A*(1-A)*(phi+N_pos))/((phi+1)*N_pos));
    double significance = A/std;
    statpars[0] = std;
    statpars[1] = significance;
    //fprintf(stderr,"STD %f \t SIGNIFICANCE %f \n",std,significance);
    double alpha;
    double beta;
    double Dx_var_numerator;
    double Dx_var_deumerator;
    double Dx_std;
    double Dx_std_norm;
    
    for(int i = 0; i < NUMROWS;i++){
      double Dx = A * pow((1 - q), fabs(XCOL[i])) + c;
      statpars[i+2] = Dx;
      //fprintf(stderr,"i: %d val: %f\n",i,Dx);
      // gzprintf(file,"%f \t",Dx);
    }

    for(int i = 0; i < NUMROWS;i++){
        double Dx = A * pow((1 - q), fabs(XCOL[i])) + c;
        alpha = Dx * phi;
        beta = (1 - Dx) * phi;
        Dx_var_numerator = NCOL[i]*alpha*beta*(alpha+beta+NCOL[i]);
        Dx_var_deumerator = pow((alpha+beta),2)*(alpha+beta+1);
        Dx_std = std::sqrt((Dx_var_numerator/Dx_var_deumerator));
        Dx_std_norm = Dx_std / NCOL[i];
	statpars[NUMROWS+2+i] = Dx_std_norm;
	//	fprintf(stderr,"at: %d val:%f\n",NUMROWS+2+i, Dx_std_norm);
	//        gzprintf(file,"%f",Dx_std_norm);
	
    }

}


#ifdef __WITH_MAIN__

int main(int argc,char **argv){
  const char *fname = "MycoBactBamSEOutSortMDSortN.mismatches.txt.gz";
  int nopt = 10;
  if(argc>1){
    for(int i=1;i<argc;i+=2)
      if(!strcmp(argv[i],"-file"))
	fname = strdup(argv[i+1]);
      else if(!strcmp(argv[i],"-nopt")){
	nopt= atoi(argv[i+1]);
      }else{
	fprintf(stderr,"\t-> Unknown option: %s\n",argv[i]);
	return 0;
      }
   
  }
  fprintf(stderr,"\t-> -fname: %s -nopt: %d\n",fname,nopt);
  double **dat = read1_ugly_matrix(fname);
  double pars[6] = {0.1,0.1,0.01,1000};//last one will contain the llh,and the ncall for the objective function
  optimoptim(pars,dat,nopt);
  double stats[2+2*(int)dat[0][0]];
  getstat(dat,pars,stats);

  //printit
  fprintf(stderr,"(A,q,c,phi,llh,ncall,SigmaD,Zfit)\n");
  for(int i=0;i<6;i++)
    fprintf(stderr,"%f\t",pars[i]);
  fprintf(stderr,"%f\t%f\n",stats[0],stats[1]);

  fprintf(stderr,"direction,cycle,k,n,dx,dx_conf\n");
  int nrows = (int) dat[0][0];
  int ncycle = nrows/2;
  double *dx = stats+2;
  double *dx_conf = stats+2+nrows;
  for(int i=0;i<ncycle;i++)
    fprintf(stderr,"5\t%d\t%.0f\t%0.f\t%f\t%f\n",i,dat[1][i],dat[2][i],dx[i],dx_conf[i]);
  dx = stats+2+ncycle;
  dx_conf = stats+2+nrows+ncycle;
  for(int i=0;i<ncycle;i++)
    fprintf(stderr,"3\t%d\t%.0f\t%0.f\t%f\t%f\n",i,dat[1][i+ncycle],dat[2][i+ncycle],dx[i],dx_conf[i]);
}

#endif
