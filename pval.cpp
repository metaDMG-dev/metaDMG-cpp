#ifndef PVAL_h
#define PVAL_h
#include <math.h>
#include <gsl/gsl_cdf.h>

#define PI 3.14159265358979323846264338327950288

// Define the t-distribution PDF
double t_pdf(double x, int dof) {
    double numerator = tgamma((dof + 1) / 2.0);
    double denominator = sqrt(dof * PI) * tgamma(dof / 2.0);
    return numerator / denominator * pow(1 + (x * x) / dof, -(dof + 1) / 2.0);
}

double t_cdf(double x, int dof, double step = 0.001) {
    double result = 0.0;
    double current_x = -20.0;  // Starting from a very negative value
    while (current_x <= x) {
        result += t_pdf(current_x, dof) * step;
        current_x += step;
    }
    return result;
}

double calculate_p_value(double test_statistic, int dof) {
    double p_value = 1.0 - t_cdf(test_statistic, dof);
    return p_value;
}

double calculate_two_tailed_p_value(double t_statistic, int degrees_of_freedom, int num_segments) {
    double cdf_positive = t_cdf(fabs(t_statistic), degrees_of_freedom, num_segments);
    double cdf_negative = t_cdf(-fabs(t_statistic), degrees_of_freedom, num_segments);

    double p_value = 2.0 * fmin(cdf_positive, cdf_negative); // Two-tailed p-value

    return p_value;
}

double calculate_one_tailed_p_value(double t_statistic, int degrees_of_freedom, int num_segments) {
    double cdf_positive = t_cdf(t_statistic, degrees_of_freedom, num_segments);
    double one_tailed_p_value = 1.0 - cdf_positive; // One-tailed p-value (upper tail)

    return one_tailed_p_value;
}

double calculate_z_score(double confidence_level) {
    // Calculate the (1 - alpha/2) quantile for a two-tailed distribution
    double alpha = 1.0 - confidence_level;
    double quantile = 1.0 - alpha / 2.0;

    // Estimate the Z-score using the quantile function
    double z_score = gsl_cdf_ugaussian_Pinv(quantile);

    return z_score;
}

#endif