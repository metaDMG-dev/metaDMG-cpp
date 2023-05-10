#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

void mu_phi_to_alpha_beta(double mu, double phi, double* alpha, double* beta);

void get_priors(double* priors);

double log_beta_func(double alpha, double beta);

double log_beta(double x, double alpha, double beta);

double log_exponential(double x, double loc, double scale);

double compute_log_likelihood(const double DMGparam[], const void *dats);

double compute_log_prior(double* Combined);

double compute_log_posterior(double* CombinedParam);

void get_stat_val(double* stats,const double Combined[]);

#endif /* LIKELIHOOD_H */