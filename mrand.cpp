#include <stdlib.h>
#include <random>
#include <climits>
#include "mrand.h"

mrand_t *mrand_alloc(int type_a,long int seedval){
  mrand_t *ret =  new mrand_t();
  ret->type = type_a;

#ifdef __APPLE__
  if(ret->type==0){
    fprintf(stderr,"\t-> Problem with drand48 reentrant, will default to erand48\n");
    ret->type = 3;
  }
#else
  if(ret->type==0)
    srand48_r(seedval,(struct drand48_data *) &ret->buf0);
#endif
  if(ret->type==1){
    ret->eng = std::default_random_engine(seedval);
    ret->distr = std::uniform_real_distribution<double>(0, 1);
    ret->distrInt = std::uniform_int_distribution<>(0,INT_MAX);
  }
  if(ret->type==2)
    ret->rand_r_seed = (unsigned int) seedval;
  if(ret->type==3){
    ret->rand_r_seed = (unsigned int) seedval;
    ret->xsubi[2] = ret->rand_r_seed >> 16;
    ret->xsubi[1] = ret->rand_r_seed & 0xffffl;
    ret->xsubi[0] = 0x330e;
  }
  if(ret->type==4){
    unsigned long long tmp = -1;
    ret->nr_inv_rec = 1/((double) tmp);//(2^64-1)^-1
    ret->nr_uvw[1] = 4101842887655102017LL;
    ret->nr_uvw[2] = 1LL;
    ret->nr_uvw[0] = seedval ^ ret->nr_uvw[1];ret->nr_int64();
    ret->nr_uvw[1] =     ret->nr_uvw[0];ret->nr_int64();
    ret->nr_uvw[2] =     ret->nr_uvw[1];ret->nr_int64();

  }
  return ret;
}

void mrand_destroy(mrand_t *mr) {
    delete mr;
}

double mrand_pop(mrand_t *mr){
  double res = -1;
  if(mr->type==0){
    #if defined(__linux__) || defined(__unix__)
    drand48_r((struct drand48_data*)&mr->buf0,&res);
#else
    fprintf(stderr,"\t-> Error: mr->type==0 is not supported on this platform, will exit\n");
    exit(1);
    #endif
  }
  else if(mr->type==1){
    res =  mr->distr(mr->eng);
  }
  else if(mr->type==2){
    int randr = rand_r(&mr->rand_r_seed);
    if (randr == RAND_MAX)
      res = (double) (randr-1)/RAND_MAX;     
    else
      res = (double) randr/RAND_MAX;
  }
  else if(mr->type==3){
    res = erand48(mr->xsubi);
  }
  else if(mr->type==4){
    res = mr->nr_inv_rec * mr->nr_int64();
  }
  else{
    fprintf(stderr,"Random parameter %d is not supported\n",mr->type);
    exit(1);
  }
  if (!(res >= 0.0 && res < 1.0)) {
    fprintf(stderr, "mrand_pop: value out of range: %.17g (type=%d)\n", res, mr->type);
    exit(1);
  }
  
  return res;
}


#ifdef __WITH_MAIN__

int main(int argc,char **argv){
  mrand_t *myrand;
  int type =0;
  int seed =1;
  long long nitems =10;
  double dogeom = 0.0;
  if(argc>3){
    type = atoi(argv[1]);
    seed = atoi(argv[2]);
    nitems = atoll(argv[3]);
  }
  if(argc>4)
    dogeom = atof(argv[4]);
  fprintf(stderr,"type: %d seed: %d nitems: %lld dogeom: %f\n",type,seed,nitems,dogeom);
  myrand = mrand_alloc(type,seed);

  if(dogeom==0){
    double sum = 0;
    for(long long i=0;i<nitems;i++){
      sum += mrand_pop(myrand);
    }
    fprintf(stdout,"type %d\tseed: %d\tsum:%f\tnitems in mio:%f mean: %f\n",type,seed,sum,nitems/1e6,sum/((double)nitems));
  }
  
  return 0;
}

#endif

//g++ mrand.cpp -std=c++11 -lm -lz -D__WITH_MAIN__ -o Rand
