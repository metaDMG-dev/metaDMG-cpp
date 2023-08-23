#ifndef PVAL_h
#define PVAL_h


// Define the t-distribution PDF
double t_pdf(double x, int dof);

// Calculate the t-distribution CDF using the trapezoidal rule
double t_cdf(double x, int dof, int num_segments);

double calculate_p_value(double test_statistic, int dof);

double calculate_two_tailed_p_value(double t_statistic, int degrees_of_freedom, int num_segments);

double calculate_one_tailed_p_value(double t_statistic, int degrees_of_freedom, int num_segments);

double calculate_z_score(double confidence_level);

#endif