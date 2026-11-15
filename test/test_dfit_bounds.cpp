// Unit test: the dfit optimiser must stay in the valid parameter region A + c < 1.
//
// D(x) = A(1-q)^x + c is largest at x = 0, where D(0) = A + c. The beta-binomial
// likelihood needs beta = (1 - D) * phi > 0, i.e. A + c < 1. The box bounds alone
// (A < 1, c <= 0.25) allow A + c up to 1.25, and on data where every C at position 0
// reads as T on top of a clear background, the fit used to end at A = 1, c = 0.106,
// with D(0) = 1.106 and a positive log-likelihood (+106.6), which is impossible.
//
// It also checks the likelihood-ratio test against the no-damage model (optim_null):
// lrt = 2 * (llh - llh_null) must be >= 0, small without damage and large but finite
// when damage is close to 1 (where Zfit = A / sigma explodes because sigma -> 0).

#include <cmath>
#include <cstdio>
#include <vector>

#include "../dfit_optim.h"
#include "../mrand.h"

static const int HOWMANY = 15;
static const int NOPT = 10;

// Build the dat layout used by main_dfit.cpp: dat[0] = {nrows, ncall, x...}, dat[1] = k, dat[2] = N.
static double **make_dat(const std::vector<double> &k, const std::vector<double> &n) {
    int nrows = 2 * HOWMANY;
    double **dat = new double *[4];
    dat[0] = new double[nrows + 2];
    dat[1] = new double[nrows];
    dat[2] = new double[nrows];
    dat[3] = new double[nrows];
    dat[0][0] = nrows;
    dat[0][1] = 0;
    for (int end = 0; end < 2; ++end)  // same counts at the 5' and the 3' end
        for (int i = 0; i < HOWMANY; ++i) {
            int row = end * HOWMANY + i;
            dat[0][row + 2] = i;
            dat[1][row] = k[i];
            dat[2][row] = n[i];
            dat[3][row] = n[i] > 0 ? k[i] / n[i] : 0.0;
        }
    return dat;
}

static void free_dat(double **dat) {
    for (int i = 0; i < 4; ++i) delete[] dat[i];
    delete[] dat;
}

static int check(const char *name, bool ok, const char *what, double got) {
    if (ok) return 0;
    std::fprintf(stderr, "[%s] %s (got %g)\n", name, what, got);
    return 1;
}

static int run_case(const char *name, const std::vector<double> &k, const std::vector<double> &n,
                    double A_exp, double q_exp, double c_exp, double tol_A, double tol_q, double tol_c,
                    double lrt_min, double lrt_max) {
    double **dat = make_dat(k, n);
    mrand_t *rng = mrand_alloc(3, 1);
    double pars[6] = {0.1, 0.1, 0.01, 1000, 0, 0};  // same start values as main_dfit.cpp
    optimoptim(pars, dat, NOPT, rng);
    double nullpars[2];
    double llh_null = optim_null(dat, nullpars);
    if (pars[4] < llh_null) {
        double fallback[6] = {};
        optim_from_null_start(fallback, dat, nullpars);
        if (fallback[4] > pars[4])
            std::memcpy(pars, fallback, sizeof(pars));
    }

    int nrows = 2 * HOWMANY;
    std::vector<double> stats(4 + 2 * nrows);
    getstat(dat, pars, stats.data());

    double A = pars[0], q = pars[1], c = pars[2], llh = pars[4];
    double raw_lrt = 2 * (llh - llh_null);
    double lrt = raw_lrt > 0.0 ? raw_lrt : 0.0;
    std::printf("[%s] A=%.4f q=%.4f c=%.4f phi=%.1f llh=%.4f A+c=%.6f llh_null=%.4f lrt=%.3f\n", name, A, q, c,
                pars[3], llh, A + c, llh_null, lrt);

    int fail = 0;
    fail += check(name, A + c < 1.0, "A + c must be < 1", A + c);
    fail += check(name, llh <= 0.0, "log-likelihood of count data must be <= 0", llh);
    fail += check(name, std::fabs(A - A_exp) < tol_A, "A far from expected value", A);
    if (q_exp >= 0)
        fail += check(name, std::fabs(q - q_exp) < tol_q, "q far from expected value", q);
    fail += check(name, std::fabs(c - c_exp) < tol_c, "c far from expected value", c);
    fail += check(name, llh_null <= 0.0, "null log-likelihood must be <= 0", llh_null);
    fail += check(name, lrt >= lrt_min && lrt <= lrt_max, "lrt outside expected range", lrt);
    for (int i = 0; i < nrows; ++i) {
        double dx = stats[4 + i], dx_conf = stats[4 + nrows + i];
        fail += check(name, dx >= 0.0 && dx < 1.0, "fitted D(x) must be in [0,1)", dx);
        fail += check(name, !std::isnan(dx_conf), "dx_conf must not be NaN", dx_conf);
    }

    mrand_destroy(rng);
    free_dat(dat);
    return fail;
}

int main() {
    int fail = 0;

    // Every C at position 0 is a T, plus 10% background everywhere: MLE is A = 0.9, c = 0.1.
    std::vector<double> k_full(HOWMANY, 100), n_full(HOWMANY, 1000);
    k_full[0] = 1000;
    fail += run_case("full_damage_pos0", k_full, n_full, 0.9, -1, 0.1, 0.01, 0, 0.005, 100, 1000);

    // Realistic strong damage that was always inside the valid region: A = 0.6, q = 0.4, c = 0.02.
    std::vector<double> k_real(HOWMANY), n_real(HOWMANY, 1000);
    for (int i = 0; i < HOWMANY; ++i)
        k_real[i] = std::round(1000 * (0.6 * std::pow(0.6, i) + 0.02));
    fail += run_case("realistic_damage", k_real, n_real, 0.6, 0.4, 0.02, 0.01, 0.01, 0.005, 100, 100000);

    // No damage at all, 2% background everywhere: A ~ 0 and the test must find nothing.
    std::vector<double> k_none(HOWMANY, 20), n_none(HOWMANY, 1000);
    fail += run_case("no_damage", k_none, n_none, 0.0, -1, 0.02, 0.01, 0, 0.005, 0, 1);

    if (fail) {
        std::fprintf(stderr, "dfit bounds test: %d check(s) failed\n", fail);
        return 1;
    }
    std::puts("dfit bounds test passed");
    return 0;
}
