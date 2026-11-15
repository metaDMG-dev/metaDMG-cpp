#include "mrand.h"

double optimoptim(double *invec,double **dat,int nopt,mrand_t *rand_alloc);
double optim_null(double **dat, double outpars[] = nullptr);
double optim_from_null_start(double *invec, double **dat, const double nullpars[]);
void getstat(double **dat,double *pars,double *statpars);
