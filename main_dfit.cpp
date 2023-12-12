/*
  16 july, this file contains functionality from the print_ugly function and codeparts in the misc.
  The function will implement the ML estimate of the dfit from metaDMG-cpp bioxarhive paper.
  
 */

#include <random>
#include <cstdio>
#include <htslib/hts.h>   // for htsFormat, hts_opt_add, htsFile, hts_opt
#include <htslib/sam.h>   // for htsFormat, hts_opt_add, htsFile, hts_opt
#include <htslib/bgzf.h>
#include <ctime>
#include <sys/time.h>
#include <math.h>
#include <gsl/gsl_cdf.h>
#include <algorithm> // for std::sort

#include <iostream>
#include "profile.h"
#include "shared.h"
#include "ngsLCA.h" //<- print_chain
#include "types.h"       // for int2intvec, int2int
#include "version.h"     // for METADAMAGE_VERSION
#include "dfit_optim.h"
#include "pval.h"
#include "dfit_helppage.h"
#include "mrand.h"

extern htsFormat *dingding2;

// the compare function for double values
static int compare (const void * a, const void * b)
{
  if (*(double*)a > *(double*)b) return 1;
  else if (*(double*)a < *(double*)b) return -1;
  else return 0;  
}

typedef struct{
  std::map<int, mydataD> *retmap;
  int howmany;
  sam_hdr_t *hdr;
  int2char *name_map;
  int libprep;
  int nopt;
  int nbootstrap;
  double CI; int doCI;
  int sigtype;
  int seed;
  int rng_type;
  int doboot;
  kstring_t *kstr;
  kstring_t *bootkstr;
  int showfits;
}ding;

void make_bootstrap_data(double **in,double **out,int howmany,mrand_t *rand_alloc){

  double *xcol_in = in[0];
  double *kvec_in = in[1];
  double *nvec_in = in[2];
  double *xcol_out = out[0];
  double *kvec_out = out[1];
  double *nvec_out = out[2];

  xcol_out[0] = xcol_in[0];
  xcol_out[1] = xcol_in[1];
  for(int p = 0;p < 2*howmany;p++){
    xcol_out[p+2] = xcol_in[p+2];
    double prob =(double) kvec_in[p]/nvec_in[p];
    //    fprintf(stderr,"Prob: %f\n",prob);
    kvec_out[p] = nvec_out[p] = 0;
    for(int i =0;i<(int)nvec_in[p];i++){
      nvec_out[p] +=1;
      if(mrand_pop(rand_alloc)<prob)
	      kvec_out[p] +=1;
    }
    //   fprintf(stderr,"k: %f n: %f\n",kvec_out[p],nvec_out[p]);
  }
}

//aa,ac,ag,at,ca,cc,cg,ct,ga
//ct and ga has index 7,8 when zero indexed
void make_dfit_format(mydataD &md,double **dat,int howmany,int libprep){
  /*
  for double-stranded lib prep we would observe both C>T and G>A
  for single-stranded lib prep we would only observe C>T

  Therefore the columns extracted from the MisMatchMatrix will differ
  */

  int col_5p = 7; // ct for ds
  int col_3p = 8; // ga for da
  if(libprep == 1){
    col_3p = 7; //only ct for ss
  }

  dat[0][0] = 2*howmany;
  dat[0][1] = 0;
  for(int i=0;i<howmany;i++){
    dat[0][i+2] =i;//plug in position

    dat[2][i] = 0;//initialize kcol to zero, a few lines down we will sum over all C*
    dat[1][i] = md.fwD[i*16+col_5p];//plugin ct at kcol
    
    for(int at=0;at<4;at++)
      dat[2][i] = dat[2][i]+md.fwD[i*16+4+at];

    dat[3][i] = (double) dat[1][i]/dat[2][i];

  }
  // In the 3' end, for ds we use the GA col (8), for ss we use CA col (7)
  for(int i=0;i<howmany;i++){
    dat[0][howmany+i+2] =i;//plug in position

    dat[2][howmany+i] = 0;//initialize kcol to zero, a few lines dows we will sum over all G*
    dat[1][howmany+i] = md.bwD[i*16+col_3p];//plugin ga at kcol
    
    for(int at=0;at<4;at++)
      dat[2][howmany+i] = dat[2][howmany+i]+md.bwD[i*16+col_3p+at];

    dat[3][howmany+i] = (double) dat[1][howmany+i]/dat[2][howmany+i];
  }

}



//aa,ac,ag,at,ca,cc,cg,ct,ga
//ct and ga has index 7,8 when zero indexed
/*
  |dat[0]| = 2*howmany+2; twice the cyclelength first is 5' next is 3'
  |dat[1]| = 2*howmany; kvec
  |dat[2]| = 2*howmany; nvec
  |dat[3]| = 2*howmany; freq(k) = kvec/nvec
 */
void make_dfit_format_bootstrap(mydataD &md,double **dat,int howmany,mrand_t *rand_alloc){
  dat[0][0] = 2*howmany;
  dat[0][1] = 0;

  //do 5'
  for(int i=0;i<howmany;i++){
    dat[0][i+2] =i;

    double tsum[4] = {0,0,0,0};
    for(int r=0;r<4;r++)//loop over ref
      for(int o=0;o<4;o++)//loop over obs
	      tsum[r] += md.fwD[i*16+r*4+o];
    //  fprintf(stderr,"tsum: %f %f\n",tsum[1],tsum[2]);
    double prob1 = tsum[1]/(tsum[0]+tsum[1]+tsum[2]+tsum[3]);//prob of ref=c
    double prob2 =((double) md.fwD[i*16+7])/tsum[1];//prob of obs=t, with ref=c
    //  fprintf(stderr,"probl1: %f %f\n",prob1,prob2);
    dat[2][i] = 0;//initialize kcol to zero, a few lines down we will sum over all C*

    for(int f=0;f<tsum[0]+tsum[1]+tsum[2]+tsum[3];f++){
      if(mrand_pop(rand_alloc)<prob1){//if reference is c
        if(mrand_pop(rand_alloc)<prob2)//if observertion is t
          dat[1][i] += 1;
	
	      dat[2][i] += 1;//we are at ref=c, so increment the nvec
      }
    }
    
    dat[3][i] = (double) dat[1][i]/dat[2][i];
  }

  //do 3'
  for(int i=0;i<howmany;i++){
    dat[0][howmany+i+2] =i;

    double tsum[4] = {0,0,0,0};
    for(int r=0;r<4;r++)
      for(int o=0;o<4;o++)
	      tsum[r] += md.bwD[i*16+r*4+o];
    
    double prob1 = tsum[2]/(tsum[0]+tsum[1]+tsum[2]+tsum[3]);//sum over all
    double prob2 =((double) md.bwD[i*16+8])/tsum[2];//ga over sum of g*
    dat[2][howmany+i] = 0;//initialize kcol to zero, a few lines dows we will sum over all G*
    for(int f=0;f<tsum[0]+tsum[1]+tsum[2]+tsum[3];f++){
      if(mrand_pop(rand_alloc)<prob1){//if reference is g
        if(mrand_pop(rand_alloc)<prob2)//if observertion is a
          dat[1][howmany+i] += 1;
        
        dat[2][howmany+i] += 1;//we are at ref=g, so increment the nvec
      }
    }

    dat[3][howmany+i] = (double) dat[1][howmany+i]/dat[2][howmany+i];
  }

}

void make_dfit_format_bootstrap2(mydataD &md,double **dat,int howmany,int seed){
  static std::random_device rd;
  static std::mt19937 gen(rd());
  
  srand48(seed);
  dat[0][0] = 2*howmany;
  dat[0][1] = 0;

  //do 5'
  for(int i=0;i<howmany;i++){
    dat[0][i+2] =i;

    double tsum[4] = {0,0,0,0};
    for(int r=0;r<4;r++)//loop over ref
      for(int o=0;o<4;o++)//loop over obs
	      tsum[r] += md.fwD[i*16+r*4+o];
    //  fprintf(stderr,"tsum: %f %f\n",tsum[1],tsum[2]);
    double prob1 = tsum[1]/(tsum[0]+tsum[1]+tsum[2]+tsum[3]);//prob of ref=c
    double prob2 =((double) md.fwD[i*16+7])/tsum[1];//prob of obs=t, with ref=c
    // fprintf(stderr,"probl1: %f %f\n",prob1,prob2);
    std::binomial_distribution<> d(tsum[0]+tsum[1]+tsum[2]+tsum[3], prob1);
    dat[2][i] = d(gen);
    std::binomial_distribution<> d2(dat[2][i], prob2);
    dat[1][i] = d2(gen);
    
    dat[3][i] = (double) dat[1][i]/dat[2][i];
  }

  //do 3'
  for(int i=0;i<howmany;i++){
    dat[0][howmany+i+2] =i;

    double tsum[4] = {0,0,0,0};
    for(int r=0;r<4;r++)
      for(int o=0;o<4;o++)
	      tsum[r] += md.bwD[i*16+r*4+o];
    
    double prob1 = tsum[2]/(tsum[0]+tsum[1]+tsum[2]+tsum[3]);//sum over all
    double prob2 =((double) md.bwD[i*16+8])/tsum[2];//ga over sum of g*

    std::binomial_distribution<> d(tsum[0]+tsum[1]+tsum[2]+tsum[3], prob1);
    
    dat[2][howmany+i] = d(gen);
    std::binomial_distribution<> d2(dat[2][i], prob2);
    dat[2][howmany+i] = d2(gen);


    dat[3][howmany+i] = (double) dat[1][howmany+i]/dat[2][howmany+i];
  }

}

mydataD getval_full(std::map<int, mydataD> &retmap, int2intvec &child, int taxid, int howmany);
mydata2 getval_stats(std::map<int, mydata2> &retmap, int2intvec &child, int taxid) ;


void make_dfit_header(kstring_t *kstr,int showfits,int nbootstrap,int howmany ){
 
  if(showfits==0){
    //fprintf(stderr,"INSIDE THE FIRST SHOWFITS loop 0 \n");
    // Without bootstrap
    if(nbootstrap < 2){
      ksprintf(kstr,"taxid\tA\tq\tc\tphi\tllh\tncall\tsigmaD\tZfit\n");
    }
    else{
      ksprintf(kstr,"taxid\tA\tq\tc\tphi\tllh\tncall\tsigmaD\tZfit\tA_b\tq_b\tc_b\tphi_b\tA_CI_l\tA_CI_h\tq_CI_l\tq_CI_h\tc_CI_l\tc_CI_h\tphi_CI_l\tphi_CI_h\n");
    }
  }
  else if(showfits==1){
    //fprintf(stderr,"INSIDE THE FIRST SHOWFITS loop 1 \n");
    // With bootstrap
      
    if(nbootstrap < 2){
      ksprintf(kstr,"taxid\tA\tq\tc\tphi\tllh\tncall\tsigmaD\tZfit");
    }
    else{
      ksprintf(kstr,"taxid\tA\tq\tc\tphi\tllh\tncall\tsigmaD\tZfit\tA_b\tq_b\tc_b\tphi_b\tA_CI_l\tA_CI_h\tq_CI_l\tq_CI_h\tc_CI_l\tc_CI_h\tphi_CI_l\tphi_CI_h");
    }
    // And fwd + bwd dx and Conf information
    for(int i=0;i<howmany;i++){
      //fprintf(stderr,"\t asfafs %d\n",i);
	    ksprintf(kstr,"\tfwdx%d\tfwdxConf%d",i,i);
    }
    for(int i=0;i<howmany;i++)
	    ksprintf(kstr,"\tbwdx%d\tbwdxConf%d",i,i);
    
    ksprintf(kstr,"\n");
  }
  else{
    //fprintf(stderr,"INSIDE THE FIRST SHOWFITS loop 2 \n");
    // With bootstrap
    if(nbootstrap < 2){
      ksprintf(kstr,"taxid\tA\tq\tc\tphi\tllh\tncall\tsigmaD\tZfit");
    }
    else{
      ksprintf(kstr,"taxid\tA\tq\tc\tphi\tllh\tncall\tsigmaD\tZfit\tA_b\tq_b\tc_b\tphi_b\tA_CI_l\tA_CI_h\tq_CI_l\tq_CI_h\tc_CI_l\tc_CI_h\tphi_CI_l\tphi_CI_h");
    }
    // And fwd + bwd k, N, dx, f and Conf information
    for(int i=0;i<howmany;i++){
      //fprintf(stderr,"\t asfafs %d\n",i);
	    ksprintf(kstr,"\tfwK%d\tfwN%d\tfwdx%d\tfwf%d\tfwdxConf%d",i,i,i,i,i);
    }
    for(int i=0;i<howmany;i++)
	    ksprintf(kstr,"\tbwK%d\tbwN%d\tbwdx%d\tbwf%d\tbwdxConf%d",i,i,i,i,i);
    ksprintf(kstr,"\n");
  }
}

void slave_block(std::map<int, mydataD> &retmap,int howmany,sam_hdr_t *hdr,int2char &name_map,int libprep,int nopt,int nbootstrap,double CI, int doCI,int sigtype,int seed,int rng_type,int doboot,kstring_t *kstr,kstring_t *bootkstr,int showfits);

void *slaveslave(void *ptr){
  ding *dng = (ding *)ptr;
  slave_block(*(dng->retmap),dng->howmany,dng->hdr,*dng->name_map,dng->libprep,dng->nopt,dng->nbootstrap,dng->CI,dng->doCI,dng->sigtype,dng->seed,dng->rng_type,dng->doboot,dng->kstr,dng->bootkstr,dng->showfits);//<-fill in the rest from the struct
  return NULL;
}

void slave_block(std::map<int, mydataD> &retmap,int howmany,sam_hdr_t *hdr,int2char &name_map,int libprep,int nopt,int nbootstrap,double CI, int doCI,int sigtype,int seed,int rng_type,int doboot,kstring_t *kstr,kstring_t *bootkstr,int showfits){
  mrand_t *rand_alloc = mrand_alloc(rng_type,seed);
  int npars = 4;

  double **dat = new double*[npars];
  dat[0] = new double[2+2*howmany];
  dat[1] = new double [2*howmany];
  dat[2] = new double [2*howmany];
  dat[3] = new double [2*howmany];

  for (std::map<int, mydataD>::iterator it = retmap.begin(); it != retmap.end(); it++) {
    int taxid = it->first;
    mydataD md = it->second;
    if (it->second.nreads == 0)
      continue;
    
    if(hdr!=NULL){
      ksprintf(kstr, "%s\t", sam_hdr_tid2name(hdr, it->first));
    }
    else{
      ksprintf(kstr,"%d\t",it->first);
    }
    
    make_dfit_format(md,dat,howmany,libprep);
    //fprintf(stderr,"after make_dfit_format\n");
    double pars[6] = {0.1,0.1,0.01,1000};//last one will contain the llh,and the ncall for the objective function
    optimoptim(pars,dat,nopt,rand_alloc);
      
    double pars_b[6] = {0.1,0.1,0.01,1000};//last one will contain the llh,and the ncall for the objective function
    
    double cistat[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
    double square_diff[6] = {0,0,0,0,0,0};
    double margin_error[6] = {0,0,0,0,0,0};
      
    double **CI_val = new double*[npars];
    
    if(nbootstrap>1){
	    //Bootstrapping procedure
      //do bootstrap, print to screen for now
      double **bootdata = new double*[npars];
      bootdata[0] = new double[2+2*howmany];
      bootdata[1] = new double [2*howmany];
      bootdata[2] = new double [2*howmany];
      bootdata[3] = new double [2*howmany];
      
      double z_score = calculate_z_score(CI);
      
      double** bootcidata = (double**)malloc(nbootstrap * sizeof(double*));
      
      double acumtmp = 0; double qcumtmp = 0; double ccumtmp = 0; double phicumtmp = 0;
      double astdtmp = 0; double qstdtmp = 0; double cstdtmp = 0; double phistdtmp = 0;
	
      //store the confidence interval  values
      for (int j = 0; j < npars; j++){
        CI_val[j] = new double [2];
      }
      
      // allocate
      for (int i = 0; i < nbootstrap; i++){
        bootcidata[i] = (double*)malloc(npars * sizeof(double));
        for (int j = 0; j < npars; j++) {
          bootcidata[i][j] = 0.0; // You can replace this with your initialization values
        }
      }
	
      //fprintf(stdout,"Bootstrap estimates\n");
      for(int b=0;b<nbootstrap;b++){
        if(sigtype==1){
          mrand_t *rand_alloc_boot = mrand_alloc(rng_type,seed+b);
          make_bootstrap_data(dat,bootdata,howmany,rand_alloc_boot); //(int)(seed+b)
        }
        else if(sigtype==2){
          mrand_t *rand_alloc_boot = mrand_alloc(rng_type,seed+b);
          make_dfit_format_bootstrap(md,bootdata,howmany,rand_alloc_boot);//<-this makes a new dat that has been resampled
        }
        else if(sigtype==3){
          make_dfit_format_bootstrap2(md,bootdata,howmany,(int)(seed+b));//<-this makes a new dat that has been resampled
        }
        
        pars_b[0] = 0.1;pars_b[1]=0.1,pars_b[2]=0.01,pars_b[3]=1000,pars_b[4]=0,pars_b[5]=0;
        optimoptim(pars_b,bootdata,nopt,rand_alloc);

        for(int ii=0;ii<npars;ii++){
          bootcidata[b][ii] = pars_b[ii];
        }
        
        //fprintf(stdout,"\n");    
        acumtmp += bootcidata[b][0];
        qcumtmp += bootcidata[b][1];
        ccumtmp += bootcidata[b][2];
        phicumtmp += bootcidata[b][3];
        
        //fprintf(stdout,"%f\n",pars[5]);
        if(doboot>0){
          if(hdr!=NULL){
            ksprintf(bootkstr, "%s", sam_hdr_tid2name(hdr, it->first));
          }
          else{
            if(name_map.size()==0)
        ksprintf(bootkstr,"%d",it->first);
            else{
        int2char::iterator nit = name_map.find(it->first);
        if(nit==name_map.end()){
          fprintf(stderr,"\t-> Problem finding taxid: %ds\n",it->first);
          exit(1);
        }
        ksprintf(bootkstr,"%d:%s",it->first,nit->second);
            }
          }
          for (int ii = 0; ii < npars; ii++) {
            //fprintf(stdout, "%d \t %f\t", b,bootcidata[b][ii]);
            ksprintf(bootkstr,"\t%f",bootcidata[b][ii]);
          }
          //ksprintf(bootkstr,"\t%f",bootcidata[b][0]-bootcidata[b][2]); //A-c
          ksprintf(bootkstr,"\n");
        }
      }
      for(int i=0;i<npars;i++){
        delete[] bootdata[i];
      }
      delete[] bootdata;  // Free the array of pointers*/

      cistat[0] = acumtmp/nbootstrap;
      cistat[1] = qcumtmp/nbootstrap;
      cistat[2] = ccumtmp/nbootstrap;
      cistat[3] = phicumtmp/nbootstrap;
	
      if(doCI == 1){
        for (int j = 0; j < npars; j++) {
          for (int i = 0; i < nbootstrap; i++) {
            square_diff[i] += (bootcidata[i][j] - cistat[j])*(bootcidata[i][j] - cistat[j]);
          }
          cistat[j+npars] = sqrt(square_diff[j] / (nbootstrap - 1));
          margin_error[j] = z_score * (cistat[j+npars]);
        }
        
        for (int j = 0; j < npars; j++){
          CI_val[j][0] = cistat[j]-margin_error[j];
          CI_val[j][1] = cistat[j]+margin_error[j];
          //std::cout << CI_val[j][0] << " sad " << CI_val[j][1] << std::endl;
        }
        
      }
      else if(doCI == 2){
        //for (int i = 0; i < nbootstrap; i++) {std::cout << bootcidata[i][0] << std::endl;}
        // Sort each column (index j) within rows independently
        for (int j = 0; j < npars; j++) {
          double* temp_column = (double*)malloc(nbootstrap * sizeof(double));
          for (int i = 0; i < nbootstrap; i++) {
            temp_column[i] = bootcidata[i][j];
          }
          qsort(temp_column, nbootstrap, sizeof(double), compare);
          for (int i = 0; i < nbootstrap; i++) {
            bootcidata[i][j] = temp_column[i];
          }
          
          free(temp_column);
        }
        //std::cout << "-----" << std::endl;
        //for (int i = 0; i < nbootstrap; i++){std::cout << bootcidata[i][0] << std::endl;}
        int CI_cutval = ceil(nbootstrap * ((1-CI)/2));
        
        for (int j = 0; j < npars; j++){
          CI_val[j][0] = bootcidata[0+CI_cutval][j]; 
          CI_val[j][1] = bootcidata[nbootstrap-CI_cutval-1][j];
          //std::cout << CI_val[j][0] << " sad " << CI_val[j][1] << std::endl;
        }
      }
      
      //for(int i=0;i<4;i++){delete[] CI_val[i];}
	
    	// Free allocated memory
      for (int i = 0; i < nbootstrap; i++){free(bootcidata[i]);}
      free(bootcidata);
    }
    
    // Sigma and Z
    double stats[2+2*(int)dat[0][0]];
    getstat(dat,pars,stats);
      
    // stats contains standard deviation, then significance, then the calculated Dx for each position then the normalized 
    ksprintf(kstr,"%f",pars[0]);
    for(int i=1;i<6;i++){
	    ksprintf(kstr,"\t%f",pars[i]);
    }
      
    //std::cout << "stats " << stats[0] << " " << stats[1] << std::endl;
    ksprintf(kstr,"\t%f\t%f",stats[0],stats[1]);
      
    if(showfits == 0){
      if(nbootstrap > 1){
        //fprintf(stderr,"INSIDE nboot loop \n");
        // adding the estimated parameters for the bootstrapping method -> this does not conform with the original estimated A value obtained from 
        // the first optimoptit(pars,dat,nopt), so this is the mean estimate across all bootstrap, same for the other parameters q,c,phi,llh,ncall
        for(int i=0;i<6;i++){
          if(i == 0){
            // A
            ksprintf(kstr,"\t%f",cistat[i]);
          }
          else if(i == 4){
            // llh
            //ksprintf(kstr,"%f\t",(-1)*cistat[i]);
            continue;
          }
          else if(i == 5){
            // ncall
            //ksprintf(kstr,"%f\t",round(cistat[i]));
            continue;
          }
          else{
            // q c phi 
            ksprintf(kstr,"\t%f",cistat[i]);
          }
        }
        ksprintf(kstr,"\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",CI_val[0][0],CI_val[0][1],CI_val[1][0],CI_val[1][1],CI_val[2][0],CI_val[2][1],CI_val[3][0],CI_val[3][1]);
      }
    }
    else if(showfits == 1){
      if(nbootstrap > 1){
        for(int i=0;i<6;i++){
          if(i == 0){
            ksprintf(kstr,"\t%f",cistat[i]);
          }
          else if(i == 4){continue;}
          else if(i == 5){continue;}
          else{
            ksprintf(kstr,"\t%f",cistat[i]);
          }
        }
        ksprintf(kstr,"\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",CI_val[0][0],CI_val[0][1],CI_val[1][0],CI_val[1][1],CI_val[2][0],CI_val[2][1],CI_val[3][0],CI_val[3][1]);
      }
      //fprintf(stderr,"-----------------\n");
      int nrows = (int) dat[0][0];
      int ncycle = nrows/2;
      double *dx = stats+2;
      double *dx_conf = stats+2+nrows;
      //std::cout << stats[3] << std::endl;
      // beginning positions fwK0	fwN0	fwdx0	fwdxConf0
      for(int i=0;i<ncycle;i++){
        if(isnan(dx_conf[i])){dx_conf[i] = 0.0;}
        if(isnan(dat[3][i])){dat[3][i] = 0.0;}
        ksprintf(kstr,"\t%f\t%f",dx[i],dx_conf[i]);
      }
      
      dx = stats+2+ncycle;
      dx_conf = stats+2+nrows+ncycle;
      
      for(int i=0;i<ncycle;i++){
        if (isnan(dx_conf[i])){dx_conf[i] = 0.0;}
        ksprintf(kstr,"\t%f\t%f",dx[i],dx_conf[i]);
      }
    }
    else if(showfits == 2){
      if(nbootstrap > 1){
        for(int i=0;i<6;i++){
          if(i == 0){
            ksprintf(kstr,"\t%f",cistat[i]);
          }
          else if(i == 4){continue;}
          else if(i == 5){continue;}
          else{
            ksprintf(kstr,"\t%f",cistat[i]);
          }
        }
        ksprintf(kstr,"\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",CI_val[0][0],CI_val[0][1],CI_val[1][0],CI_val[1][1],CI_val[2][0],CI_val[2][1],CI_val[3][0],CI_val[3][1]);
      }

      //fprintf(stderr,"-----------------\n");
      int nrows = (int) dat[0][0];
      int ncycle = nrows/2;
      double *dx = stats+2;
      double *dx_conf = stats+2+nrows;
      //std::cout << stats[2] << " " << dx[0] << std::endl;
          
      // beginning positions fwK0	fwN0	fwdx0	fwdxConf0
      for(int i=0;i<ncycle;i++){
        //std::cout << " check " << dat[1][i] << " " << dat[2][i]<< " " <<dat[3][i]<< " " <<dx[i]<< " " <<dx_conf[i] << std::endl;
        if(isnan(dx_conf[i])){dx_conf[i] = 0.0;}
        if(isnan(dat[3][i])){dat[3][i] = 0.0;}
        ksprintf(kstr,"\t%.0f\t%0.f\t%f\t%f\t%f",dat[1][i],dat[2][i],dat[3][i],dx[i],dx_conf[i]);
      }
	
      dx = stats+2+ncycle;
      dx_conf = stats+2+nrows+ncycle;
      
      for(int i=0;i<ncycle;i++){
        if (isnan(dx_conf[i])){dx_conf[i] = 0.0;}
        if(isnan(dat[3][i+ncycle])){dat[3][i+ncycle] = 0.0;} //the f column
        ksprintf(kstr,"\t%.0f\t%0.f\t%f\t%f\t%f",dat[1][i+ncycle],dat[2][i+ncycle],dat[3][i+ncycle],dx[i],dx_conf[i]);
      }
    }
    ksprintf(kstr,"\n");
    delete[] CI_val;  // Free the array of pointers

  }

  // Free the memory
  for(int i=0;i<npars;i++){
    delete[] dat[i];
  }
  delete[] dat;  // Free the array of pointers
}



int main_dfit(int argc, char **argv) {
    /*
    fprintf(stderr, "./metaDMG-cpp dfit file.bdamage.gz -names file.gz -nodes trestructure.gz -lcastat fil.gz -bam file.bam -showfits int -nopt int -nbootstrap int -seed int -doCI int -CI float -lib <ds,ss> -out file\n");
    fprintf(stderr, "-------------\n Estimate damage patterns with beta-binomial model\n");
    fprintf(stderr, "\tEstimate damage patterns for each chr/scaffold contig (local mode), using lca stats\n");
    fprintf(stderr, "\t\t./metaDMG-cpp dfit file.bdamage.gz -names file.gz -nodes trestructure.gz -lcastat fil.gz -bam file.bam -showfits int -nopt int -nbootstrap int -seed int -doCI int -CI float -lib <ds,ss> -out file\n");
    fprintf(stderr, "\tEstimate one global damage pattern \n");
    fprintf(stderr, "-------------\n Estimate damage patterns with binomial model\n");
    fprintf(stderr, "\tEstimate damage patterns for each chr/scaffold contig (local mode), using lca stats\n");
    fprintf(stderr, "\t\t./metaDMG-cpp dfit file.bdamage.gz -names file.gz -nodes trestructure.gz -lcastat fil.gz -bam file.bam -showfits int -nopt int -nbootstrap int -seed int -doCI int -CI float -lib <ds,ss> -out file\n");
    fprintf(stderr, "\tEstimate one global damage pattern \n");
    fprintf(stderr, "\t\t./metaDMG-cpp dfit metaDMG-cpp/metaDMG-cpp dfit Pitch6getDMG.bdamage.gz -doboot 1 -nbootstrap 5 -nopt 10 -showfits 0\n");
    */
    if (argc <= 1){
      HelpPageSimple(stderr);
      return 0;
    }
    char *infile_bdamage = NULL;
    char *infile_nodes = NULL;
    char *infile_names = NULL;
    char *infile_bam = NULL;
    char *outfile_name = NULL;
    char *lib_prep = NULL;
    int howmany;//this is the cycle
    int showfits=0;
    int nopt = 10;
    int sigtype = 1;
    int nbootstrap = 1;
    int doboot = 0;
    int seed = time(NULL);
    int nthreads = 1;
    double CI = 0.95;
    int doCI = 2;
    int rng_type = -1;
    while (*(++argv)) {
        if (strcasecmp("-h", *argv) == 0)
          HelpPageSimple(stderr);
        else if (strcasecmp("--help", *argv) == 0)
          HelpPage(stderr);
        else if (strcasecmp("--names", *argv) == 0)
            infile_names = strdup(*(++argv));
        else if (strcasecmp("--nodes", *argv) == 0)
            infile_nodes = strdup(*(++argv));
        else if (strcasecmp("--bam", *argv) == 0)
            infile_bam = strdup(*(++argv));
        else if (strcasecmp("--out", *argv) == 0 || strcasecmp("--out_prefix", *argv) == 0)
            outfile_name = strdup(*(++argv));
        else if (strcasecmp("--nopt", *argv) == 0)
          nopt = atoi(*(++argv));
        else if (strcasecmp("--sigtype", *argv) == 0)
          sigtype = atoi(*(++argv));
        else if (strcasecmp("--showfits", *argv) == 0)
          showfits = atoi(*(++argv));
        else if (strcasecmp("--nbootstrap", *argv) == 0)
          nbootstrap = atoi(*(++argv));
        else if (strcasecmp("--doboot", *argv) == 0)
          doboot = atoi(*(++argv));
        else if (strcasecmp("--seed", *argv) == 0)
          seed = atoi(*(++argv));
        else if(strcasecmp("-rng",*argv)==0 || strcasecmp("--rand",*argv)==0){
          rng_type = atoi(*(++argv));
        }
        else if (strcasecmp("--CI", *argv) == 0)
          CI = atof(*(++argv));
        else if (strcasecmp("--doCI", *argv) == 0)
          doCI = atoi(*(++argv));
        else if (strcasecmp("--lib", *argv) == 0)
            lib_prep = strdup(*(++argv));
	      else if (strcasecmp("--nthreads", *argv) == 0)
            nthreads = atoi(*(++argv));
        else
          infile_bdamage = strdup(*argv);
    }
    if(infile_nodes&&!infile_names){
      fprintf(stderr,"\t-> -names file.gz must be defined with -nodes is defined\n");
      exit(1);
    }

    // pseudo random number generator specific to OS.
    if (rng_type == -1){
      #if defined(__linux__) || defined(__unix__)
        rng_type = 0;
      #elif defined(__APPLE__) || defined(__MACH__)
        rng_type = 3;
        //when 0 it will have problems with drand48 reentrant, will default to erand48 (MacroRandType 3)
      #else
      #   error "Unknown compiler"
      #endif
    }


    clock_t t = clock();
    struct timeval start_time, end_time;
    gettimeofday(&start_time, NULL);

    htsFile *samfp = NULL;
    sam_hdr_t *hdr = NULL;
    if (infile_bam) {
        if ((samfp = sam_open_format(infile_bam, "r", dingding2)) == NULL) {
            fprintf(stderr, "[%s] nonexistant file: %s\n", __FUNCTION__, infile_bam);
            exit(1);
        }
        hdr = sam_hdr_read(samfp);
    }
    if(infile_bam!=NULL&&infile_nodes!=NULL){
      fprintf(stderr,"\t-> Potential issue, supplying both -bam and -nodes might not be meaningfull\n");
      exit(1);
    }
    
    //fprintf(stderr, "infile_names: %s infile_bdamage: %s nodes: %s lca_stat: %s infile_bam: %s showfits: %d nopt: %d outname: %s sigtype: %d nbootstrap: %d", infile_names, infile_bdamage, infile_nodes, infile_lcastat, infile_bam,showfits,nopt,outfile_name,sigtype,nbootstrap);
    //fprintf(stderr, "#VERSION:%s\n", METADAMAGE_VERSION);
    if(outfile_name==NULL)
      outfile_name = strdup(infile_bdamage);
    char buf[1024];
    snprintf(buf, 1024, "%s.dfit.gz", outfile_name);
    fprintf(stderr, "\t-> Dumping file: \'%s\'\n", buf);
    BGZF *fpfpfp = bgzf_open(buf, "wb");
   

    char bootbuf[1024];
    BGZF *bootfp;
    kstring_t *bootkstr = new kstring_t;
    bootkstr->s = NULL; bootkstr->l = bootkstr->m = 0;
    if(doboot>0){
      snprintf(bootbuf, 1024, "%s.boot.stat.gz", outfile_name);
      bootfp = bgzf_open(bootbuf, "wb");    
      ksprintf(bootkstr,"taxid\tA_b\tq_b\tc_b\tphi_b\n");
    }

    // map of taxid -> taxid
    int2int parent;
    // map of taxid -> rank
    int2char rank;
    // map of parent -> child taxids
    int2intvec child;

    if (infile_nodes != NULL)
        parse_nodes(infile_nodes, rank, parent, child, 1);

    std::map<int, mydataD> retmap = load_bdamage_full(infile_bdamage, howmany);
    fprintf(stderr, "\t-> Number of entries in damage pattern file: %lu printlength(howmany):%d\n", retmap.size(), howmany);

    int2char name_map;

    if (infile_names)
        name_map = parse_names(infile_names);

    float presize = retmap.size();
    if(child.size()>0)
      getval_full(retmap, child, 1, howmany);  // this will do everything
    float postsize = retmap.size();
    fprintf(stderr, "\t-> pre: %f post:%f grownbyfactor: %f\n", presize, postsize, postsize / presize);

    int libprep = 0;

    if (lib_prep != NULL){
      if(strcasecmp(lib_prep, "ds")==0){
        libprep = 0;
      }
      else if(strcasecmp(lib_prep, "ss")==0){
        libprep = 1;
      }
    }

    //prepare matrix to be passed to optimization
    kstring_t *kstr = new kstring_t;
    kstr->s = NULL; kstr->l = kstr->m = 0;
    fprintf(stderr,"\t-> Will do optimization of %lu different taxids/chromosomes/scaffolds\n",retmap.size());
    make_dfit_header(kstr,showfits,nbootstrap,howmany);


    {//loop over threads, for now we have no threads
      if(nthreads==1){
      kstring_t *kstr_block = new kstring_t;
      kstr_block->s = NULL; kstr_block->l = kstr_block->m = 0;
      
      kstring_t *bootkstr_block = new kstring_t;
      bootkstr_block->s = NULL; bootkstr_block->l = bootkstr_block->m = 0;
      slave_block(retmap,howmany,hdr,name_map,libprep,nopt,nbootstrap,CI,doCI,sigtype,seed,rng_type,doboot,kstr_block,bootkstr_block,showfits);
      
      ksprintf(kstr,"%s",kstr_block->s);
      ksprintf(bootkstr,"%s",bootkstr_block->s);

      free(kstr_block->s);
      delete kstr_block;
      free(bootkstr_block->s);
      delete bootkstr_block;
      }
      else{
        std::map<int, mydataD> *ary = new std::map<int, mydataD>[nthreads];
        int cnt = 0;
        for (std::map<int, mydataD>::iterator it = retmap.begin(); it != retmap.end(); it++){
          ary[cnt % nthreads] [it->first] = it->second;
          cnt++;
        }
	
	      ding *dings = new ding[nthreads];
        for(int i=0;i<nthreads;i++){
          dings[i].retmap = &(ary[i]);
          dings[i].howmany = howmany;
          dings[i].hdr = hdr;
          dings[i].name_map = &name_map;
          dings[i].libprep = libprep;
          dings[i].nopt = nopt;
          dings[i].nbootstrap = nbootstrap;
          dings[i].CI = CI;
          dings[i].doCI = doCI;
          dings[i].sigtype = sigtype;
          dings[i].seed = (seed+i)*100;
          dings[i].rng_type = rng_type;
          dings[i].doboot = doboot;
          fprintf(stderr, "\t-> Initiating thread %d with thread specific seed %d inferred from global seed %d with pseudo-random number generator type  %d\n",i,dings[i].seed/100,seed,rng_type);
          //fprintf(stderr,"INITIATED THREAD WHAT %d WITH SEED VALUE WHAT %d WITH SPECIFIC SEED %d and seedtype %d \n",i,seed,dings[i].seed,rng_type);

          kstring_t *kstr = new kstring_t;
          kstr->s = NULL; kstr->l = kstr->m = 0;
          dings[i].kstr = kstr;
          kstring_t *bootkstr = new kstring_t;
          bootkstr->s = NULL; bootkstr->l = bootkstr->m = 0;
          dings[i].bootkstr = bootkstr;
          dings[i].showfits = showfits;
        }
        pthread_t mythreads[nthreads];
        for(int i=0;i<nthreads;i++)
          pthread_create(&mythreads[i],NULL,slaveslave, &dings[i]);
        
        for(int i=0;i<nthreads;i++)
          pthread_join(mythreads[i],NULL);
          
        for(int i=0;i<nthreads;i++){
          ksprintf(kstr,"%s",dings[i].kstr->s);
          ksprintf(bootkstr,"%s",dings[i].bootkstr->s);
        }
      }
    }
    
    if(bgzf_write(fpfpfp,kstr->s,kstr->l) == 0){
      fprintf(stderr, "\t-> Cannot write to output BGZ file\n");
      exit(1);
    }
    kstr->l = 0;
    bgzf_close(fpfpfp);
    
    
    if(doboot>0){
      if(bgzf_write(bootfp,bootkstr->s,bootkstr->l) == 0){
        fprintf(stderr, "\t-> Cannot write to output BGZ file\n");
        exit(1);
      }
      bootkstr->l = 0;
      bgzf_close(bootfp);
    }
    
    for(int2char::iterator it=name_map.begin();it!=name_map.end();it++)
      free(it->second);
    for(int2char::iterator it=rank.begin();it!=rank.end();it++)
      free(it->second);
    
    for( std::map<int, mydataD>::iterator it = retmap.begin();it!=retmap.end();it++){
      mydataD md = it->second;
      delete [] md.fwD;
      delete [] md.bwD;
    }
    
    free(kstr->s);
    delete kstr;

    free(bootkstr->s);
    delete bootkstr;
    
    if (hdr)
        bam_hdr_destroy(hdr);
    if (samfp)
        sam_close(samfp);
    if(infile_bdamage)
      free(infile_bdamage);
    if(infile_nodes)
      free(infile_nodes);
    if(infile_names)
      free(infile_names);
    if(infile_bam)
      free(infile_bam);
    if(outfile_name)
      free(outfile_name);
    if(lib_prep)
      free(lib_prep);
    
  gettimeofday(&end_time, NULL);
  long seconds = end_time.tv_sec - start_time.tv_sec;
  long microseconds = end_time.tv_usec - start_time.tv_usec;
  double elapsed_time = seconds + microseconds / 1000000.0;

  fprintf(stderr, "\t[ALL done] cpu-time used =  %.6f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.6f sec\n", elapsed_time);
  return 0;
}
