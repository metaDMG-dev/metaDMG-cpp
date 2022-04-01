
#include <__utility/move.h>         // for move
#include <gsl/gsl_vector_double.h>  // for gsl_vector

#include <__config>  // for std
#include <cstdlib>   // for size_t
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <fstream>  // for string
#include <new>      // for operator new
#include <vector>   // for vector
//#include "matplotlibcpp.h"
using namespace std;
using namespace Eigen;
// namespace plt = matplotlibcpp;

struct par {
    MatrixXd X;
    MatrixXd Tab;
    int order;
    int str_pt;
    int end_pt;
    int numpos;
};

struct res {
    int dir;
    int refnucid;
    double LR_like;
    double MLR_like;
    string refnuc;
    VectorXd LR_coeff;
    VectorXd MLR_coeff;
    MatrixXd X;
    vector<vector<double> > LR_freq_p;
    vector<vector<double> > LR_freq_n;
    vector<vector<double> > MLR_freq_p;
    vector<vector<double> > MLR_freq_n;
};

typedef struct {
    // filenames
    const char *namedir;
    const char *outname;
    // const char* outfigname;
    const char *outfreqname;

    int model;
    int numppos;
    int numnpos;
    int order;
} pars;

void readdata(const char *filename, string *ColumnName, size_t **Table);
void matrixselector(size_t **Table, int numppos, int numnpos, int str_pt, int end_pt, int dir, VectorXd &b, MatrixXd &Tab, double &scale);
void constructX(MatrixXd &X, int numpos, int numcol, int order);
double loglike(const gsl_vector *v, void *params);
double calloglike(VectorXd &R_coeff, MatrixXd &X, MatrixXd &Tab, int order, int str_pt, int end_pt, int numpos);
int gslminloglike(MatrixXd &X, MatrixXd &Tab, int order, int str_pt, int end_pt, int numpos, VectorXd LR_coeff, VectorXd &MLR_coeff);
void modelcalculation(size_t **Table, int numppos, int numnpos, int model, int order, vector<res> &Results, string &figname);
void Freqcalculator4cond(vector<res> &Results, int numppos, int numnpos);
void Freqcalculator4uncon(vector<res> &Results, int numppos, int numnpos);
// void freqplt(vector<res > &Results, size_t ** Table, int model, int numppos, int numnpos, int order, double like, string &figname);
void output(const char *outputname, const char *filename, string figname, vector<res> &Results, int model, int order, double like);
pars *get_pars(int argc, char **argv);
pars *pars_init();
int main_regression(int argc, char **argv);
