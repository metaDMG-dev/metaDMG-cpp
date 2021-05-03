#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <gsl/gsl_types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include "matplotlibcpp.h"
#include "regression.h"
using namespace std;
using namespace Eigen;
namespace plt = matplotlibcpp;

//Get the table
void readdata(const char* filename, string* ColumnName,size_t** Table){
    ifstream infile;
    infile.open(filename);
    string line;
    if(infile.is_open()){
        getline(infile, line);
        istringstream iss(line);
        string word;
        int i = 0;
        while(iss >> word)
        {
            ColumnName[i++] = word;
        }
        int row = 0;
        size_t tab;
        while(getline(infile, line)){
            int col = 0;
            istringstream iss1(line);
            while (getline(iss1,word,' ')){
                int n = word.size();
                if (n>0){
                    if (word[n-1]!='\''){
                        stringstream ss2(word);
                        ss2 >> tab;
                        Table[row][col] = tab;
                    }else{
                        stringstream ss3(word.substr(0,n-1));
                        ss3 >> tab;
                        Table[row][col] = tab;
                    }
                    col = col + 1;
                }
            }
            row = row + 1;
        }
        
    }
    infile.close();
}

void matrixselector(size_t** Table, int numppos, int numnpos, int str_pt, int end_pt, int dir, VectorXd &b, MatrixXd &Tab, double &scale){
    int k = end_pt - str_pt;
    scale = (double)Table[0][4];
    cout << "Scale is "<<scale<<"\n";
    if (dir > 0){
        for (int i = 0; i < numppos; i++){
            int nj = 0;
            double sum = 0;
            for (int j = str_pt; j< end_pt; j++){
                b(nj+k*i) = log((double)Table[i][j]/(double)Table[i][end_pt]);
                //cout << b(nj+k*i) << "\t";
                Tab(i,nj)=(double)Table[i][j]/scale;
                sum += Tab(i,nj);
                nj = nj + 1;
            }
            Tab(i,k) = sum+(double)Table[i][end_pt]/scale;
            //cout << "\n";
            //cout << Table[i][end_pt] << "\n";
        }
    }else if(dir < 0){
        int ni = 0;
        for (int i = numppos; i < numppos+numnpos; i++){
            int nj = 0;
            double sum = 0;
            for (int j = str_pt; j< end_pt; j++){
                //cout << Table[i][j] << "\t";
                b(nj+k*ni) = log((double)Table[i][j]/(double)Table[i][end_pt]);
                Tab(ni,nj)=(double)Table[i][j]/scale;
                sum += Tab(ni,nj);
                //cout << b(nj+k*ni) << "\t";
                nj = nj + 1;
            }
            Tab(ni,k) = sum+(double)Table[i][end_pt]/scale;
            ni = ni + 1;
            //cout<<"\n";
        }
    }else{
        for (int i = 0; i < numppos; i++){
            int nj = 0;
            double sum = 0;
            for (int j = str_pt; j< end_pt; j++){
                //cout << Table[i][j] << "\t";
                b(nj+k*i) = log((double)Table[i][j]/(double)Table[i][end_pt]);
                Tab(i,nj)=(double)Table[i][j]/scale;
                sum += Tab(i,nj);
                //                cout << b(nj+k*i) << "\n";
                cout << Tab(i,nj) << "\t";
                nj = nj + 1;
            }
            Tab(i,k) = sum+(double)Table[i][end_pt]/scale;
            //cout<<"\n";
            cout << Tab(i,k) << "\n";
        }
        for (int i = numppos; i < numppos+numnpos; i++){
            int nj = 0;
            double sum = 0;
            for (int j = str_pt; j< end_pt; j++){
                //cout << Table[i][23-j] << "\t";
                b(nj+k*i) = log((double)Table[i][23-j]/(double)Table[i][23-end_pt]);
                Tab(i,nj)=(double)Table[i][23-j]/scale;
                sum += Tab(i,nj);
                //cout << Table[i][j] << "\t";
                cout << Tab(i,nj) << "\t";
                nj = nj + 1;
            }
            Tab(i,k) = sum+(double)Table[i][23-end_pt]/scale;
            cout << Tab(i,k) << "\n";
        }
    }
}

void constructX(MatrixXd &X, int numpos, int numcol, int order){
    //    for (int l = 0; l <= order; l++){
    //        for(int i=0; i<numpos; i++){
    //            for(int j=0; j < numcol; j++){
    //                X(i*numcol+j,l*numcol+j)=(i+1)^l;??
    //            }
    //        }
    //    }
    for (int l=0; l<numcol; l++){
        for(int i=0; i<numpos; i++){
            for(int j=0; j<=order; j++){
                X(l+numcol*i,l+j*numcol)=pow((i+1),j);
            }
        }
    }
    
}

//struct par{
//    MatrixXd X;
//    MatrixXd Tab;
//    int order;
//    int str_pt;
//    int end_pt;
//    int numpos;
//};
//
//struct res{
//    int dir;
//    int refnucid;
//    double LR_like;
//    double MLR_like;
//    string refnuc;
//    VectorXd LR_coeff;
//    VectorXd MLR_coeff;
//    MatrixXd X;
//    vector<vector<double> > LR_freq_p;
//    vector<vector<double> > LR_freq_n;
//    vector<vector<double> > MLR_freq_p;
//    vector<vector<double> > MLR_freq_n;
//};
// Likelihood function
//double loglike(const gsl_vector *v, void *params)
//{
//    double l = 0;
//    par *p = (par *)params;
//    MatrixXd X = p->X;
//    MatrixXd Tab = p->Tab;
//    int typenum = p->end_pt-p->str_pt;
//    int len = typenum*(p->order+1);
//    int numpos = p->numpos;
//    VectorXd vec(len);
//
//    for (int i=0; i < len; i++){
//        vec(i) = gsl_vector_get(v, i);
//    }
//    VectorXd Y = X*vec;
//    double ymax = Y.maxCoeff();
//    for (int i = 0; i < numpos; i++){
//        double sumE = 0;
//        for (int j = 0; j < typenum; j++){
//            l = l + Tab(i,j)*Y(i*typenum+j);
//            sumE = sumE + exp(Y(i*typenum+j)-ymax);
//        }
//        sumE = sumE + exp(-ymax);
//        l = l + Tab(i,typenum)*(-ymax-log(sumE));
//    }
//    return -l;
//}

double loglike(const gsl_vector *v, void *params)
{
    par *p = (par *)params;
    MatrixXd X = p->X;
    MatrixXd Tab = p->Tab;
    int typenum = p->end_pt-p->str_pt;
    int len = typenum*(p->order+1);
    int numpos = p->numpos;
    int Y_len = typenum*numpos;
    VectorXd vec(len);
    
    for (int i=0; i < len; i++){
        vec(i) = gsl_vector_get(v, i);
    }
    
    double l = 0;
    VectorXd Y = X*vec;
    for (int k=0; k<Y_len; k=k+typenum){
        double sumEBX = 0;
        double sumBX = 0;
        int pos = k/typenum;
        for (int j=0;j<typenum;j++){
            sumEBX = sumEBX + exp(Y(k+j));
            sumBX = sumBX + Tab(pos,j)*Y(k+j);
        }
        l = l - Tab(pos,typenum)*log(1+sumEBX)+sumBX;
    }
    return -l;
}

double calloglike(VectorXd &R_coeff, MatrixXd &X, MatrixXd &Tab, int order, int str_pt, int end_pt, int numpos)
{
    int typenum = end_pt-str_pt;
    int len = typenum*(order+1);
    int Y_len = typenum*numpos;
    VectorXd vec(len);
    for (int i=0; i < len; i++){
        vec(i) = R_coeff(i);
    }
    double l = 0;
    VectorXd Y = X*vec;
    for (int k=0; k<Y_len; k=k+typenum){
        double sumEBX = 0;
        double sumBX = 0;
        int pos = k/typenum;
        for (int j=0;j<typenum;j++){
            sumEBX = sumEBX + exp(Y(k+j));
            sumBX = sumBX + Tab(pos,j)*Y(k+j);
        }
        l = l - Tab(pos,typenum)*log(1+sumEBX) + sumBX;
    }
    return l;
}

int gslminloglike (MatrixXd &X, MatrixXd &Tab, int order, int str_pt, int end_pt, int numpos, VectorXd LR_coeff, VectorXd &MLR_coeff)
{
    //double par[5] = {1.0, 2.0, 10.0, 20.0, 30.0};
    par* params = new par;
    params->X = X;
    params->order = order;
    params->str_pt = str_pt;
    params->end_pt = end_pt;
    params->numpos = numpos;
    params->Tab = Tab;
    int typenum = end_pt-str_pt;
    int len = typenum*(order+1);
    
    //Check
    //cout<<"likelihood of LR_coeff2 is " <<loglike2(LR_coeff, params)<<"\n";
    
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minloglike;
    
    size_t iter = 0;
    int status;
    double size;
    
    /* Starting point */
    x = gsl_vector_alloc(len);
    for (int i = 0; i < len; i++){
        gsl_vector_set (x, i, LR_coeff(i));
    }
    
    /* Set initial step sizes to 1 */
    ss = gsl_vector_alloc(len);
    gsl_vector_set_all (ss, 1e-5);
    
    /* Initialize method and iterate */
    minloglike.n = len;
    minloglike.f = loglike;
    minloglike.params = params;
    
    s = gsl_multimin_fminimizer_alloc(T, len);
    gsl_multimin_fminimizer_set(s, &minloglike, x, ss);
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if (status)
            break;
        
        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-6);
        
        if (status == GSL_SUCCESS)
        {
            printf ("converged to minimum at\n");
            printf ("%lu %10.3e %10.3e f() = %7.3f size = %.3f\n",
                    iter,
                    gsl_vector_get (s->x, 0),
                    gsl_vector_get (s->x, 1),
                    s->fval, size);
        }
        
    }
    while (status == GSL_CONTINUE && iter < 10000);
    for (int i = 0; i < len; i++){
        MLR_coeff(i)=gsl_vector_get(s->x, i);
    }
    
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    
    //Check
    //cout<<"likelihood of MLR_coeff is " <<loglike1(MLR_coeff, params)<<"\n";
    //cout<<"likelihood of MLR_coeff2 is " <<loglike2(MLR_coeff, params)<<"\n";
    return status;
}

void modelcalculation(size_t** Table, int numppos, int numnpos, int model, int order, vector<res > &Results, string &figname){
    if (model<=3){
        string NucName = "ACGT";
        double scale;
        double like;
        if (model == 0){
            cout << "Unfolded full regression is under conducted!\n";
            figname = "Unfolded full regression";
            int str_pt = 4;
            int end_pt = 19;
            int dir = 1;
            int b_len = 15*numppos;
            VectorXd bp(b_len);
            MatrixXd Tabp = MatrixXd::Zero(numppos,16);
            cout << "Dir is " << dir <<"\n";
            matrixselector(Table, numppos, numnpos, str_pt, end_pt,dir,bp,Tabp,scale);
            MatrixXd X = MatrixXd::Zero(b_len,15*(order+1));
            constructX(X, numppos, 15, order);
            VectorXd LR_coeff_p = X.bdcSvd(ComputeThinU | ComputeThinV).solve(bp);
            VectorXd MLR_coeff_p = VectorXd::Zero(15*(order+1));
            gslminloglike(X, Tabp, order, str_pt, end_pt, numppos, LR_coeff_p, MLR_coeff_p);
            res tmp;
            tmp.dir = dir;
            tmp.refnucid = 4;
            tmp.refnuc = "N";
            tmp.LR_coeff = LR_coeff_p;
            tmp.MLR_coeff = MLR_coeff_p;
            tmp.X = X;
            like = calloglike(LR_coeff_p, X, Tabp, order, str_pt, end_pt, numppos);
            tmp.LR_like = scale*like;
            like = calloglike(MLR_coeff_p, X, Tabp, order, str_pt, end_pt, numppos);
            tmp.MLR_like = scale*like;
            Results.push_back(tmp);
            
            
            dir = -1;
            b_len = 15*numnpos;
            VectorXd bn(b_len);
            MatrixXd Tabn = MatrixXd::Zero(numnpos,16);
            cout << "Dir is " << dir <<"\n";
            matrixselector(Table, numppos, numnpos, str_pt, end_pt,dir,bn, Tabn,scale);
            MatrixXd Y = MatrixXd::Zero(15*numnpos,15*(order+1));
            constructX(Y, numppos, 15, order);
            VectorXd LR_coeff_n = Y.bdcSvd(ComputeThinU | ComputeThinV).solve(bn);
            VectorXd MLR_coeff_n = VectorXd::Zero(15*(order+1));
            gslminloglike(X, Tabn, order, str_pt, end_pt, numppos, LR_coeff_n, MLR_coeff_n);
            tmp.dir = dir;
            tmp.refnucid = 4;
            tmp.refnuc = "N";
            tmp.LR_coeff = LR_coeff_n;
            tmp.MLR_coeff = MLR_coeff_n;
            tmp.X = Y;
            like = calloglike(LR_coeff_n, Y, Tabn, order, str_pt, end_pt, numnpos);
            tmp.LR_like = scale*like;
            like = calloglike(MLR_coeff_n, Y, Tabn, order, str_pt, end_pt, numnpos);
            tmp.MLR_like = scale*like;
            Results.push_back(tmp);
        }else if(model == 1){
            cout << "Unfolded conditional regression is under conducted!\n";
            figname = "Unfolded conditional regression";
            for (int i=0; i<4; i++){
                int str_pt = 4*i+4;
                int end_pt = 4*i+7;
                
                int dir = 1;
                int b_len = 3*numppos;
                cout << "Reference is " << NucName[i] << ", dir is " << dir <<"\n";
                VectorXd bp(b_len);
                MatrixXd Tabp = MatrixXd::Zero(numppos,4);
                matrixselector(Table, numppos, numnpos, str_pt, end_pt,dir,bp,Tabp,scale);
                MatrixXd X = MatrixXd::Zero(b_len,3*(order+1));
                constructX(X, numppos, 3, order);
                VectorXd LR_coeff_p = X.bdcSvd(ComputeThinU | ComputeThinV).solve(bp);
                VectorXd MLR_coeff_p = VectorXd::Zero(3*(order+1));
                cout<<Tabp<<"\n";
                gslminloglike(X, Tabp, order, str_pt, end_pt, numppos, LR_coeff_p, MLR_coeff_p);
                res tmp;
                tmp.dir = dir;
                tmp.refnucid = i;
                tmp.refnuc = NucName[i];
                tmp.LR_coeff = LR_coeff_p;
                tmp.MLR_coeff = MLR_coeff_p;
                tmp.X = X;
                like = calloglike(LR_coeff_p, X, Tabp, order, str_pt, end_pt, numppos);
                tmp.LR_like = scale*like;
                like = calloglike(MLR_coeff_p, X,Tabp, order, str_pt, end_pt, numppos);
                tmp.MLR_like = scale*like;
                Results.push_back(tmp);
                
                dir = -1;
                b_len = 3*numnpos;
                cout << "Reference is " << NucName[i] << ", dir is " << dir <<"\n";
                VectorXd bn(b_len);
                MatrixXd Tabn = MatrixXd::Zero(numnpos,4);
                matrixselector(Table, numppos, numnpos, str_pt, end_pt, dir, bn, Tabn,scale);
                MatrixXd Y = MatrixXd::Zero(b_len,3*(order+1));
                constructX(Y, numnpos, 3, order);
                VectorXd LR_coeff_n = Y.bdcSvd(ComputeThinU | ComputeThinV).solve(bn);
                VectorXd MLR_coeff_n = VectorXd::Zero(3*(order+1));
                cout<<Tabn<<"\n";
                gslminloglike(X, Tabn, order, str_pt, end_pt, numnpos, LR_coeff_n, MLR_coeff_n);
                tmp.dir = dir;
                tmp.refnucid = i;
                tmp.refnuc = NucName[i];
                tmp.LR_coeff = LR_coeff_n;
                tmp.MLR_coeff = MLR_coeff_n;
                tmp.X = Y;
                like = calloglike(LR_coeff_n, Y, Tabn, order, str_pt, end_pt, numnpos);
                tmp.LR_like = scale*like;
                like = calloglike(MLR_coeff_n,Y, Tabn, order, str_pt, end_pt, numnpos);
                tmp.MLR_like = scale*like;
                cout<<"Test likelihood is "<<scale<<"\n";
                Results.push_back(tmp);
                
                cout<<"Results size "<<Results.size()<<"\n";
            }
        }else if(model == 2){
            cout << "Folded full regression is under conducted!\n";
            figname = "Folded full regression";
            int dir = 0;
            cout << "Dir is " << dir <<"\n";
            int str_pt = 4;
            int end_pt = 19;
            int b_len = 15*(numppos+numnpos);
            VectorXd b(b_len);
            MatrixXd Tab = MatrixXd::Zero(numppos+numnpos,16);
            matrixselector(Table, numppos, numnpos, str_pt, end_pt, dir, b, Tab,scale);
            MatrixXd X = MatrixXd::Zero(15*numppos,15*(order+1));
            constructX(X, numppos, 15, order);
            MatrixXd Y = MatrixXd::Zero(15*numnpos,15*(order+1));
            constructX(Y, numnpos, 15, order);
            MatrixXd Z(X.rows()+Y.rows(), X.cols());
            Z << X, Y;
            VectorXd LR_coeff = Z.bdcSvd(ComputeThinU | ComputeThinV).solve(b);
            VectorXd MLR_coeff = VectorXd::Zero(15*(order+1));
            gslminloglike(X, Tab, order, str_pt, end_pt, numppos, LR_coeff, MLR_coeff);
            res tmp;
            tmp.dir = dir;
            tmp.refnucid = 4;
            tmp.refnuc = "N";
            tmp.LR_coeff = LR_coeff;
            tmp.MLR_coeff = MLR_coeff;
            tmp.X = Z;
            like = calloglike(LR_coeff, Z, Tab, order, str_pt, end_pt, numppos+numnpos);
            tmp.LR_like = scale*like;
            like = calloglike(MLR_coeff, Z, Tab, order, str_pt, end_pt, numppos+numnpos);
            tmp.MLR_like = scale*like;
            Results.push_back(tmp);
        }else{
            cout << "Folded conditional regression is under conducted!\n";
            figname = "Folded conditional regression";
            int dir = 0;
            for (int i=0; i<4; i++){
                int str_pt = 4*i+4;
                int end_pt = 4*i+7;
                int b_len = 3*(numppos+numnpos);
                VectorXd b(b_len);
                MatrixXd Tab = MatrixXd::Zero(numppos+numnpos,4);
                cout << "Reference is " << NucName[i] << ", dir is " << dir <<"\n";
                matrixselector(Table, numppos, numnpos, str_pt, end_pt, dir, b, Tab,scale);
                MatrixXd X = MatrixXd::Zero(3*numppos,3*(order+1));
                constructX(X, numppos, 3, order);
                MatrixXd Y = MatrixXd::Zero(3*numnpos,3*(order+1));
                constructX(Y, numnpos, 3, order);
                MatrixXd Z(X.rows()+Y.rows(), X.cols());
                Z << X, Y;
                VectorXd LR_coeff = Z.bdcSvd(ComputeThinU | ComputeThinV).solve(b);
                VectorXd MLR_coeff = VectorXd::Zero(3*(order+1));
                gslminloglike(Z, Tab, order, str_pt, end_pt, numppos+numnpos, LR_coeff, MLR_coeff);
                for (int i = 0; i<3*(order+1);i++)
                cout << LR_coeff(i) << MLR_coeff(i)<<"\n";
                res tmp;
                tmp.dir = dir;
                tmp.refnucid = i;
                tmp.refnuc = NucName[i];
                tmp.LR_coeff = LR_coeff;
                tmp.MLR_coeff = MLR_coeff;
                tmp.X = Z;
                like = calloglike(LR_coeff, Z, Tab, order, str_pt, end_pt, numppos+numnpos);
                tmp.LR_like = scale*like;
                like = calloglike(MLR_coeff, Z, Tab, order, str_pt, end_pt, numppos+numnpos);
                tmp.MLR_like = scale*like;
                Results.push_back(tmp);
            }
        }
    }else if(model == 4){
        cout << "Briggs regression is under conducted! Under construction now!\n";
    }else{
        cout << "Please specify a model!\n";
        exit(-1);
    }
}
//
void Freqcalculator4cond(vector<res> &Results, int numppos, int numnpos){
    for (int i = 0; i < Results.size(); i++){
        int dir = Results[i].dir;
        int refnucid = Results[i].refnucid;
        MatrixXd X = Results[i].X;
        VectorXd LR_coeff = Results[i].LR_coeff;
        VectorXd MLR_coeff = Results[i].MLR_coeff;
        if (dir>0){
            vector<double> freq0(numppos);
            vector<double> freq1(numppos);
            vector<double> freq2(numppos);
            vector<double> freq3(numppos);
            int typenum = round(X.rows()/numppos);
            VectorXd LR_b = X*LR_coeff;
            for (int k=0; k<numppos; k++){
                freq0.at(numppos-k-1) = exp(LR_b(k*typenum));
                freq1.at(numppos-k-1) = exp(LR_b(k*typenum+1));
                freq2.at(numppos-k-1) = exp(LR_b(k*typenum+2));
                double sum = freq0.at(numppos-k-1)+freq1.at(numppos-k-1)+freq2.at(numppos-k-1)+1;
                freq0.at(numppos-k-1) = freq0.at(numppos-k-1)/sum;
                freq1.at(numppos-k-1) = freq1.at(numppos-k-1)/sum;
                freq2.at(numppos-k-1) = freq2.at(numppos-k-1)/sum;
                freq3.at(numppos-k-1) = (double)1/sum;
            }
            Results[i].LR_freq_p.push_back(freq0);
            Results[i].LR_freq_p.push_back(freq1);
            Results[i].LR_freq_p.push_back(freq2);
            Results[i].LR_freq_p.push_back(freq3);
            
            VectorXd MLR_b = X*MLR_coeff;
            for (int k=0; k<numppos; k++){
                freq0.at(numppos-k-1) = exp(MLR_b(k*typenum));
                freq1.at(numppos-k-1) = exp(MLR_b(k*typenum+1));
                freq2.at(numppos-k-1) = exp(MLR_b(k*typenum+2));
                double sum = freq0.at(numppos-k-1)+freq1.at(numppos-k-1)+freq2.at(numppos-k-1)+1;
                freq0.at(numppos-k-1) = freq0.at(numppos-k-1)/sum;
                freq1.at(numppos-k-1) = freq1.at(numppos-k-1)/sum;
                freq2.at(numppos-k-1) = freq2.at(numppos-k-1)/sum;
                freq3.at(numppos-k-1) = (double)1/sum;
            }
            Results[i].MLR_freq_p.push_back(freq0);
            Results[i].MLR_freq_p.push_back(freq1);
            Results[i].MLR_freq_p.push_back(freq2);
            Results[i].MLR_freq_p.push_back(freq3);
        }else if(dir<0){
            vector<double> freq0(numnpos);
            vector<double> freq1(numnpos);
            vector<double> freq2(numnpos);
            vector<double> freq3(numnpos);
            int typenum = round(X.rows()/numnpos);
            VectorXd LR_b = X*LR_coeff;
            for (int k=0; k<numnpos; k++){
                freq0.at(k) = exp(LR_b(k*typenum));
                freq1.at(k) = exp(LR_b(k*typenum+1));
                freq2.at(k) = exp(LR_b(k*typenum+2));
                double sum = freq0.at(k)+freq1.at(k)+freq2.at(k)+1;
                freq0.at(k) = freq0.at(k)/sum;
                freq1.at(k) = freq1.at(k)/sum;
                freq2.at(k) = freq2.at(k)/sum;
                freq3.at(k) = (double)1/sum;
            }
            Results[i].LR_freq_n.push_back(freq0);
            Results[i].LR_freq_n.push_back(freq1);
            Results[i].LR_freq_n.push_back(freq2);
            Results[i].LR_freq_n.push_back(freq3);
            
            VectorXd MLR_b = X*MLR_coeff;
            for (int k=0; k<numnpos; k++){
                freq0.at(k) = exp(MLR_b(k*typenum));
                freq1.at(k) = exp(MLR_b(k*typenum+1));
                freq2.at(k) = exp(MLR_b(k*typenum+2));
                double sum = freq0.at(k)+freq1.at(k)+freq2.at(k)+1;
                freq0.at(k) = freq0.at(k)/sum;
                freq1.at(k) = freq1.at(k)/sum;
                freq2.at(k) = freq2.at(k)/sum;
                freq3.at(k) = (double)1/sum;
            }
            Results[i].MLR_freq_n.push_back(freq0);
            Results[i].MLR_freq_n.push_back(freq1);
            Results[i].MLR_freq_n.push_back(freq2);
            Results[i].MLR_freq_n.push_back(freq3);
        }else{
            vector<double> freq0p(numppos);
            vector<double> freq1p(numppos);
            vector<double> freq2p(numppos);
            vector<double> freq3p(numppos);
            vector<double> freq0n(numnpos);
            vector<double> freq1n(numnpos);
            vector<double> freq2n(numnpos);
            vector<double> freq3n(numnpos);
            int typenum = round(X.rows()/(numnpos+numppos));
            VectorXd LR_b = X*LR_coeff;
            for (int k=0; k<numppos; k++){
                freq0p.at(numppos-k-1) = exp(LR_b(k*typenum));
                freq1p.at(numppos-k-1) = exp(LR_b(k*typenum+1));
                freq2p.at(numppos-k-1) = exp(LR_b(k*typenum+2));
                double sum = freq0p.at(numppos-k-1)+freq1p.at(numppos-k-1)+freq2p.at(numppos-k-1)+1;
                freq0p.at(numppos-k-1) = freq0p.at(numppos-k-1)/sum;
                freq1p.at(numppos-k-1) = freq1p.at(numppos-k-1)/sum;
                freq2p.at(numppos-k-1) = freq2p.at(numppos-k-1)/sum;
                freq3p.at(numppos-k-1) = (double)1/sum;
            }
            Results[i].LR_freq_p.push_back(freq0p);
            Results[i].LR_freq_p.push_back(freq1p);
            Results[i].LR_freq_p.push_back(freq2p);
            Results[i].LR_freq_p.push_back(freq3p);
            int j = 0;
            while (Results[j].refnucid != 3-refnucid){
                j++;
            }
            for (int k=numppos; k<numnpos+numnpos; k++){
                int k0 = k-numppos;
                freq0n.at(k0) = exp(LR_b(k*typenum));
                freq1n.at(k0) = exp(LR_b(k*typenum+1));
                freq2n.at(k0) = exp(LR_b(k*typenum+2));
                double sum = freq0n.at(k0)+freq1n.at(k0)+freq2n.at(k0)+1;
                freq0n.at(k0) = freq0n.at(k0)/sum;
                freq1n.at(k0) = freq1n.at(k0)/sum;
                freq2n.at(k0) = freq2n.at(k0)/sum;
                freq3n.at(k0) = (double)1/sum;
            }
            Results[j].LR_freq_n.push_back(freq3n);
            Results[j].LR_freq_n.push_back(freq2n);
            Results[j].LR_freq_n.push_back(freq1n);
            Results[j].LR_freq_n.push_back(freq0n);
            
            VectorXd MLR_b = X*MLR_coeff;
            for (int k=0; k<numppos; k++){
                freq0p.at(numppos-k-1) = exp(MLR_b(k*typenum));
                freq1p.at(numppos-k-1) = exp(MLR_b(k*typenum+1));
                freq2p.at(numppos-k-1) = exp(MLR_b(k*typenum+2));
                double sum = freq0p.at(numppos-k-1)+freq1p.at(numppos-k-1)+freq2p.at(numppos-k-1)+1;
                freq0p.at(numppos-k-1) = freq0p.at(numppos-k-1)/sum;
                freq1p.at(numppos-k-1) = freq1p.at(numppos-k-1)/sum;
                freq2p.at(numppos-k-1) = freq2p.at(numppos-k-1)/sum;
                freq3p.at(numppos-k-1) = (double)1/sum;
            }
            Results[i].MLR_freq_p.push_back(freq0p);
            Results[i].MLR_freq_p.push_back(freq1p);
            Results[i].MLR_freq_p.push_back(freq2p);
            Results[i].MLR_freq_p.push_back(freq3p);
            j = 0;
            while (Results[j].refnucid != 3-refnucid){
                j++;
            }
            for (int k=numppos; k<numnpos+numnpos; k++){
                int k0 = k-numppos;
                freq0n.at(k0) = exp(MLR_b(k*typenum));
                freq1n.at(k0) = exp(MLR_b(k*typenum+1));
                freq2n.at(k0) = exp(MLR_b(k*typenum+2));
                double sum = freq0n.at(k0)+freq1n.at(k0)+freq2n.at(k0)+1;
                freq0n.at(k0) = freq0n.at(k0)/sum;
                freq1n.at(k0) = freq1n.at(k0)/sum;
                freq2n.at(k0) = freq2n.at(k0)/sum;
                freq3n.at(k0) = (double)1/sum;
            }
            Results[j].MLR_freq_n.push_back(freq3n);
            Results[j].MLR_freq_n.push_back(freq2n);
            Results[j].MLR_freq_n.push_back(freq1n);
            Results[j].MLR_freq_n.push_back(freq0n);
        }
    }
}


void Freqcalculator4uncon(vector<res> &Results, int numppos, int numnpos){
    for (int i = 0; i < Results.size(); i++){
        int dir = Results[i].dir;
        int refnucid = Results[i].refnucid;
        MatrixXd X = Results[i].X;
        VectorXd LR_coeff = Results[i].LR_coeff;
        VectorXd MLR_coeff = Results[i].MLR_coeff;
        if (dir>0){
            vector<double> freq0(numppos);
            vector<double> freq1(numppos);
            vector<double> freq2(numppos);
            vector<double> freq3(numppos);
            vector<double> freq4(numppos);
            vector<double> freq5(numppos);
            vector<double> freq6(numppos);
            vector<double> freq7(numppos);
            vector<double> freq8(numppos);
            vector<double> freq9(numppos);
            vector<double> freq10(numppos);
            vector<double> freq11(numppos);
            vector<double> freq12(numppos);
            vector<double> freq13(numppos);
            vector<double> freq14(numppos);
            vector<double> freq15(numppos);
            int typenum = round(X.rows()/numppos);
            VectorXd LR_b = X*LR_coeff;
            for (int k=0; k<numppos; k++){
                freq0.at(numppos-k-1) = exp(LR_b(k*typenum));
                freq1.at(numppos-k-1) = exp(LR_b(k*typenum+1));
                freq2.at(numppos-k-1) = exp(LR_b(k*typenum+2));
                freq3.at(numppos-k-1) = exp(LR_b(k*typenum+3));
                freq4.at(numppos-k-1) = exp(LR_b(k*typenum+4));
                freq5.at(numppos-k-1) = exp(LR_b(k*typenum+5));
                freq6.at(numppos-k-1) = exp(LR_b(k*typenum+6));
                freq7.at(numppos-k-1) = exp(LR_b(k*typenum+7));
                freq8.at(numppos-k-1) = exp(LR_b(k*typenum+8));
                freq9.at(numppos-k-1) = exp(LR_b(k*typenum+9));
                freq10.at(numppos-k-1) = exp(LR_b(k*typenum+10));
                freq11.at(numppos-k-1) = exp(LR_b(k*typenum+11));
                freq12.at(numppos-k-1) = exp(LR_b(k*typenum+12));
                freq13.at(numppos-k-1) = exp(LR_b(k*typenum+13));
                freq14.at(numppos-k-1) = exp(LR_b(k*typenum+14));
                double sum = freq0.at(numppos-k-1)+freq1.at(numppos-k-1)+freq2.at(numppos-k-1)+freq3.at(numppos-k-1)+freq4.at(numppos-k-1)+freq5.at(numppos-k-1)+freq6.at(numppos-k-1)+freq7.at(numppos-k-1)+freq8.at(numppos-k-1)+freq9.at(numppos-k-1)+freq10.at(numppos-k-1)+freq11.at(numppos-k-1)+freq11.at(numppos-k-1)+freq12.at(numppos-k-1)+freq13.at(numppos-k-1)+freq14.at(numppos-k-1)+1;
                freq0.at(numppos-k-1) = freq0.at(numppos-k-1)/sum;
                freq1.at(numppos-k-1) = freq1.at(numppos-k-1)/sum;
                freq2.at(numppos-k-1) = freq2.at(numppos-k-1)/sum;
                freq3.at(numppos-k-1) = freq3.at(numppos-k-1)/sum;
                freq4.at(numppos-k-1) = freq4.at(numppos-k-1)/sum;
                freq5.at(numppos-k-1) = freq5.at(numppos-k-1)/sum;
                freq6.at(numppos-k-1) = freq6.at(numppos-k-1)/sum;
                freq7.at(numppos-k-1) = freq7.at(numppos-k-1)/sum;
                freq8.at(numppos-k-1) = freq8.at(numppos-k-1)/sum;
                freq9.at(numppos-k-1) = freq9.at(numppos-k-1)/sum;
                freq10.at(numppos-k-1) = freq10.at(numppos-k-1)/sum;
                freq11.at(numppos-k-1) = freq11.at(numppos-k-1)/sum;
                freq12.at(numppos-k-1) = freq12.at(numppos-k-1)/sum;
                freq13.at(numppos-k-1) = freq13.at(numppos-k-1)/sum;
                freq14.at(numppos-k-1) = freq14.at(numppos-k-1)/sum;
                freq15.at(numppos-k-1) = (double)1/sum;
            }
            Results[i].LR_freq_p.push_back(freq0);
            Results[i].LR_freq_p.push_back(freq1);
            Results[i].LR_freq_p.push_back(freq2);
            Results[i].LR_freq_p.push_back(freq3);
            Results[i].LR_freq_p.push_back(freq4);
            Results[i].LR_freq_p.push_back(freq5);
            Results[i].LR_freq_p.push_back(freq6);
            Results[i].LR_freq_p.push_back(freq7);
            Results[i].LR_freq_p.push_back(freq8);
            Results[i].LR_freq_p.push_back(freq9);
            Results[i].LR_freq_p.push_back(freq10);
            Results[i].LR_freq_p.push_back(freq11);
            Results[i].LR_freq_p.push_back(freq12);
            Results[i].LR_freq_p.push_back(freq13);
            Results[i].LR_freq_p.push_back(freq14);
            Results[i].LR_freq_p.push_back(freq15);
            
            VectorXd MLR_b = X*MLR_coeff;
            for (int k=0; k<numppos; k++){
                freq0.at(numppos-k-1) = exp(MLR_b(k*typenum));
                freq1.at(numppos-k-1) = exp(MLR_b(k*typenum+1));
                freq2.at(numppos-k-1) = exp(MLR_b(k*typenum+2));
                freq3.at(numppos-k-1) = exp(MLR_b(k*typenum+3));
                freq4.at(numppos-k-1) = exp(MLR_b(k*typenum+4));
                freq5.at(numppos-k-1) = exp(MLR_b(k*typenum+5));
                freq6.at(numppos-k-1) = exp(MLR_b(k*typenum+6));
                freq7.at(numppos-k-1) = exp(MLR_b(k*typenum+7));
                freq8.at(numppos-k-1) = exp(MLR_b(k*typenum+8));
                freq9.at(numppos-k-1) = exp(MLR_b(k*typenum+9));
                freq10.at(numppos-k-1) = exp(MLR_b(k*typenum+10));
                freq11.at(numppos-k-1) = exp(MLR_b(k*typenum+11));
                freq12.at(numppos-k-1) = exp(MLR_b(k*typenum+12));
                freq13.at(numppos-k-1) = exp(MLR_b(k*typenum+13));
                freq14.at(numppos-k-1) = exp(MLR_b(k*typenum+14));
                double sum = freq0.at(numppos-k-1)+freq1.at(numppos-k-1)+freq2.at(numppos-k-1)+freq3.at(numppos-k-1)+freq4.at(numppos-k-1)+freq5.at(numppos-k-1)+freq6.at(numppos-k-1)+freq7.at(numppos-k-1)+freq8.at(numppos-k-1)+freq9.at(numppos-k-1)+freq10.at(numppos-k-1)+freq11.at(numppos-k-1)+freq11.at(numppos-k-1)+freq12.at(numppos-k-1)+freq13.at(numppos-k-1)+freq14.at(numppos-k-1)+1;
                freq0.at(numppos-k-1) = freq0.at(numppos-k-1)/sum;
                freq1.at(numppos-k-1) = freq1.at(numppos-k-1)/sum;
                freq2.at(numppos-k-1) = freq2.at(numppos-k-1)/sum;
                freq3.at(numppos-k-1) = freq3.at(numppos-k-1)/sum;
                freq4.at(numppos-k-1) = freq4.at(numppos-k-1)/sum;
                freq5.at(numppos-k-1) = freq5.at(numppos-k-1)/sum;
                freq6.at(numppos-k-1) = freq6.at(numppos-k-1)/sum;
                freq7.at(numppos-k-1) = freq7.at(numppos-k-1)/sum;
                freq8.at(numppos-k-1) = freq8.at(numppos-k-1)/sum;
                freq9.at(numppos-k-1) = freq9.at(numppos-k-1)/sum;
                freq10.at(numppos-k-1) = freq10.at(numppos-k-1)/sum;
                freq11.at(numppos-k-1) = freq11.at(numppos-k-1)/sum;
                freq12.at(numppos-k-1) = freq12.at(numppos-k-1)/sum;
                freq13.at(numppos-k-1) = freq13.at(numppos-k-1)/sum;
                freq14.at(numppos-k-1) = freq14.at(numppos-k-1)/sum;
                freq15.at(numppos-k-1) = (double)1/sum;
            }
            Results[i].MLR_freq_p.push_back(freq0);
            Results[i].MLR_freq_p.push_back(freq1);
            Results[i].MLR_freq_p.push_back(freq2);
            Results[i].MLR_freq_p.push_back(freq3);
            Results[i].MLR_freq_p.push_back(freq4);
            Results[i].MLR_freq_p.push_back(freq5);
            Results[i].MLR_freq_p.push_back(freq6);
            Results[i].MLR_freq_p.push_back(freq7);
            Results[i].MLR_freq_p.push_back(freq8);
            Results[i].MLR_freq_p.push_back(freq9);
            Results[i].MLR_freq_p.push_back(freq10);
            Results[i].MLR_freq_p.push_back(freq11);
            Results[i].MLR_freq_p.push_back(freq12);
            Results[i].MLR_freq_p.push_back(freq13);
            Results[i].MLR_freq_p.push_back(freq14);
            Results[i].MLR_freq_p.push_back(freq15);
        }else if(dir<0){
            vector<double> freq0(numnpos);
            vector<double> freq1(numnpos);
            vector<double> freq2(numnpos);
            vector<double> freq3(numnpos);
            vector<double> freq4(numnpos);
            vector<double> freq5(numnpos);
            vector<double> freq6(numnpos);
            vector<double> freq7(numnpos);
            vector<double> freq8(numnpos);
            vector<double> freq9(numnpos);
            vector<double> freq10(numnpos);
            vector<double> freq11(numnpos);
            vector<double> freq12(numnpos);
            vector<double> freq13(numnpos);
            vector<double> freq14(numnpos);
            vector<double> freq15(numnpos);
            int typenum = round(X.rows()/numnpos);
            VectorXd LR_b = X*LR_coeff;
            for (int k=0; k<numnpos; k++){
                freq0.at(k) = exp(LR_b(k*typenum));
                freq1.at(k) = exp(LR_b(k*typenum+1));
                freq2.at(k) = exp(LR_b(k*typenum+2));
                freq3.at(k) = exp(LR_b(k*typenum+3));
                freq4.at(k) = exp(LR_b(k*typenum+4));
                freq5.at(k) = exp(LR_b(k*typenum+5));
                freq6.at(k) = exp(LR_b(k*typenum+6));
                freq7.at(k) = exp(LR_b(k*typenum+7));
                freq8.at(k) = exp(LR_b(k*typenum+8));
                freq9.at(k) = exp(LR_b(k*typenum+9));
                freq10.at(k) = exp(LR_b(k*typenum+10));
                freq11.at(k) = exp(LR_b(k*typenum+11));
                freq12.at(k) = exp(LR_b(k*typenum+12));
                freq13.at(k) = exp(LR_b(k*typenum+13));
                freq14.at(k) = exp(LR_b(k*typenum+14));
                double sum = freq0.at(k)+freq1.at(k)+freq2.at(k)+freq3.at(k)+freq4.at(k)+freq5.at(k)+freq6.at(k)+freq7.at(k)+freq8.at(k)+freq9.at(k)+freq10.at(k)+freq11.at(k)+freq11.at(k)+freq12.at(k)+freq13.at(k)+freq14.at(k)+1;
                freq0.at(k) = freq0.at(k)/sum;
                freq1.at(k) = freq1.at(k)/sum;
                freq2.at(k) = freq2.at(k)/sum;
                freq3.at(k) = freq3.at(k)/sum;
                freq4.at(k) = freq4.at(k)/sum;
                freq5.at(k) = freq5.at(k)/sum;
                freq6.at(k) = freq6.at(k)/sum;
                freq7.at(k) = freq7.at(k)/sum;
                freq8.at(k) = freq8.at(k)/sum;
                freq9.at(k) = freq9.at(k)/sum;
                freq10.at(k) = freq10.at(k)/sum;
                freq11.at(k) = freq11.at(k)/sum;
                freq12.at(k) = freq12.at(k)/sum;
                freq13.at(k) = freq13.at(k)/sum;
                freq14.at(k) = freq14.at(k)/sum;
                freq15.at(k) = (double)1/sum;
            }
            Results[i].LR_freq_n.push_back(freq0);
            Results[i].LR_freq_n.push_back(freq1);
            Results[i].LR_freq_n.push_back(freq2);
            Results[i].LR_freq_n.push_back(freq3);
            Results[i].LR_freq_n.push_back(freq4);
            Results[i].LR_freq_n.push_back(freq5);
            Results[i].LR_freq_n.push_back(freq6);
            Results[i].LR_freq_n.push_back(freq7);
            Results[i].LR_freq_n.push_back(freq8);
            Results[i].LR_freq_n.push_back(freq9);
            Results[i].LR_freq_n.push_back(freq10);
            Results[i].LR_freq_n.push_back(freq11);
            Results[i].LR_freq_n.push_back(freq12);
            Results[i].LR_freq_n.push_back(freq13);
            Results[i].LR_freq_n.push_back(freq14);
            Results[i].LR_freq_n.push_back(freq15);
            
            VectorXd MLR_b = X*MLR_coeff;
            for (int k=0; k<numnpos; k++){
                freq0.at(k) = exp(MLR_b(k*typenum));
                freq1.at(k) = exp(MLR_b(k*typenum+1));
                freq2.at(k) = exp(MLR_b(k*typenum+2));
                freq3.at(k) = exp(MLR_b(k*typenum+3));
                freq4.at(k) = exp(MLR_b(k*typenum+4));
                freq5.at(k) = exp(MLR_b(k*typenum+5));
                freq6.at(k) = exp(MLR_b(k*typenum+6));
                freq7.at(k) = exp(MLR_b(k*typenum+7));
                freq8.at(k) = exp(MLR_b(k*typenum+8));
                freq9.at(k) = exp(MLR_b(k*typenum+9));
                freq10.at(k) = exp(MLR_b(k*typenum+10));
                freq11.at(k) = exp(MLR_b(k*typenum+11));
                freq12.at(k) = exp(MLR_b(k*typenum+12));
                freq13.at(k) = exp(MLR_b(k*typenum+13));
                freq14.at(k) = exp(MLR_b(k*typenum+14));
                double sum = freq0.at(k)+freq1.at(k)+freq2.at(k)+freq3.at(k)+freq4.at(k)+freq5.at(k)+freq6.at(k)+freq7.at(k)+freq8.at(k)+freq9.at(k)+freq10.at(k)+freq11.at(k)+freq11.at(k)+freq12.at(k)+freq13.at(k)+freq14.at(k)+1;
                freq0.at(k) = freq0.at(k)/sum;
                freq1.at(k) = freq1.at(k)/sum;
                freq2.at(k) = freq2.at(k)/sum;
                freq3.at(k) = freq3.at(k)/sum;
                freq4.at(k) = freq4.at(k)/sum;
                freq5.at(k) = freq5.at(k)/sum;
                freq6.at(k) = freq6.at(k)/sum;
                freq7.at(k) = freq7.at(k)/sum;
                freq8.at(k) = freq8.at(k)/sum;
                freq9.at(k) = freq9.at(k)/sum;
                freq10.at(k) = freq10.at(k)/sum;
                freq11.at(k) = freq11.at(k)/sum;
                freq12.at(k) = freq12.at(k)/sum;
                freq13.at(k) = freq13.at(k)/sum;
                freq14.at(k) = freq14.at(k)/sum;
                freq15.at(k) = (double)1/sum;
            }
            Results[i].MLR_freq_n.push_back(freq0);
            Results[i].MLR_freq_n.push_back(freq1);
            Results[i].MLR_freq_n.push_back(freq2);
            Results[i].MLR_freq_n.push_back(freq3);
            Results[i].MLR_freq_n.push_back(freq4);
            Results[i].MLR_freq_n.push_back(freq5);
            Results[i].MLR_freq_n.push_back(freq6);
            Results[i].MLR_freq_n.push_back(freq7);
            Results[i].MLR_freq_n.push_back(freq8);
            Results[i].MLR_freq_n.push_back(freq9);
            Results[i].MLR_freq_n.push_back(freq10);
            Results[i].MLR_freq_n.push_back(freq11);
            Results[i].MLR_freq_n.push_back(freq12);
            Results[i].MLR_freq_n.push_back(freq13);
            Results[i].MLR_freq_n.push_back(freq14);
            Results[i].MLR_freq_n.push_back(freq15);
        }else{
            vector<double> freq0p(numppos);
            vector<double> freq1p(numppos);
            vector<double> freq2p(numppos);
            vector<double> freq3p(numppos);
            vector<double> freq4p(numppos);
            vector<double> freq5p(numppos);
            vector<double> freq6p(numppos);
            vector<double> freq7p(numppos);
            vector<double> freq8p(numppos);
            vector<double> freq9p(numppos);
            vector<double> freq10p(numppos);
            vector<double> freq11p(numppos);
            vector<double> freq12p(numppos);
            vector<double> freq13p(numppos);
            vector<double> freq14p(numppos);
            vector<double> freq15p(numppos);
            vector<double> freq0n(numnpos);
            vector<double> freq1n(numnpos);
            vector<double> freq2n(numnpos);
            vector<double> freq3n(numnpos);
            vector<double> freq4n(numnpos);
            vector<double> freq5n(numnpos);
            vector<double> freq6n(numnpos);
            vector<double> freq7n(numnpos);
            vector<double> freq8n(numnpos);
            vector<double> freq9n(numnpos);
            vector<double> freq10n(numnpos);
            vector<double> freq11n(numnpos);
            vector<double> freq12n(numnpos);
            vector<double> freq13n(numnpos);
            vector<double> freq14n(numnpos);
            vector<double> freq15n(numnpos);
            int typenum = round(X.rows()/(numnpos+numppos));
            VectorXd LR_b = X*LR_coeff;
            for (int k=0; k<numppos; k++){
                freq0p.at(numppos-k-1) = exp(LR_b(k*typenum));
                freq1p.at(numppos-k-1) = exp(LR_b(k*typenum+1));
                freq2p.at(numppos-k-1) = exp(LR_b(k*typenum+2));
                freq3p.at(numppos-k-1) = exp(LR_b(k*typenum+3));
                freq4p.at(numppos-k-1) = exp(LR_b(k*typenum+4));
                freq5p.at(numppos-k-1) = exp(LR_b(k*typenum+5));
                freq6p.at(numppos-k-1) = exp(LR_b(k*typenum+6));
                freq7p.at(numppos-k-1) = exp(LR_b(k*typenum+7));
                freq8p.at(numppos-k-1) = exp(LR_b(k*typenum+8));
                freq9p.at(numppos-k-1) = exp(LR_b(k*typenum+9));
                freq10p.at(numppos-k-1) = exp(LR_b(k*typenum+10));
                freq11p.at(numppos-k-1) = exp(LR_b(k*typenum+11));
                freq12p.at(numppos-k-1) = exp(LR_b(k*typenum+12));
                freq13p.at(numppos-k-1) = exp(LR_b(k*typenum+13));
                freq14p.at(numppos-k-1) = exp(LR_b(k*typenum+14));
                double sum = freq0p.at(numppos-k-1)+freq1p.at(numppos-k-1)+freq2p.at(numppos-k-1)+freq3p.at(numppos-k-1)+freq4p.at(numppos-k-1)+freq5p.at(numppos-k-1)+freq6p.at(numppos-k-1)+freq7p.at(numppos-k-1)+freq8p.at(numppos-k-1)+freq9p.at(numppos-k-1)+freq10p.at(numppos-k-1)+freq11p.at(numppos-k-1)+freq11p.at(numppos-k-1)+freq12p.at(numppos-k-1)+freq13p.at(numppos-k-1)+freq14p.at(numppos-k-1)+1;
                freq0p.at(numppos-k-1) = freq0p.at(numppos-k-1)/sum;
                freq1p.at(numppos-k-1) = freq1p.at(numppos-k-1)/sum;
                freq2p.at(numppos-k-1) = freq2p.at(numppos-k-1)/sum;
                freq3p.at(numppos-k-1) = freq3p.at(numppos-k-1)/sum;
                freq4p.at(numppos-k-1) = freq4p.at(numppos-k-1)/sum;
                freq5p.at(numppos-k-1) = freq5p.at(numppos-k-1)/sum;
                freq6p.at(numppos-k-1) = freq6p.at(numppos-k-1)/sum;
                freq7p.at(numppos-k-1) = freq7p.at(numppos-k-1)/sum;
                freq8p.at(numppos-k-1) = freq8p.at(numppos-k-1)/sum;
                freq9p.at(numppos-k-1) = freq9p.at(numppos-k-1)/sum;
                freq10p.at(numppos-k-1) = freq10p.at(numppos-k-1)/sum;
                freq11p.at(numppos-k-1) = freq11p.at(numppos-k-1)/sum;
                freq12p.at(numppos-k-1) = freq12p.at(numppos-k-1)/sum;
                freq13p.at(numppos-k-1) = freq13p.at(numppos-k-1)/sum;
                freq14p.at(numppos-k-1) = freq14p.at(numppos-k-1)/sum;
                freq15p.at(numppos-k-1) = (double)1/sum;
            }
            Results[i].LR_freq_p.push_back(freq0p);
            Results[i].LR_freq_p.push_back(freq1p);
            Results[i].LR_freq_p.push_back(freq2p);
            Results[i].LR_freq_p.push_back(freq3p);
            Results[i].LR_freq_p.push_back(freq4p);
            Results[i].LR_freq_p.push_back(freq5p);
            Results[i].LR_freq_p.push_back(freq6p);
            Results[i].LR_freq_p.push_back(freq7p);
            Results[i].LR_freq_p.push_back(freq8p);
            Results[i].LR_freq_p.push_back(freq9p);
            Results[i].LR_freq_p.push_back(freq10p);
            Results[i].LR_freq_p.push_back(freq11p);
            Results[i].LR_freq_p.push_back(freq12p);
            Results[i].LR_freq_p.push_back(freq13p);
            Results[i].LR_freq_p.push_back(freq14p);
            Results[i].LR_freq_p.push_back(freq15p);
            
            for (int k=numppos; k<numnpos+numnpos; k++){
                int k0 = k-numppos;
                freq0n.at(k0) = exp(LR_b(k*typenum));
                freq1n.at(k0) = exp(LR_b(k*typenum+1));
                freq2n.at(k0) = exp(LR_b(k*typenum+2));
                freq3n.at(k0) = exp(LR_b(k*typenum+3));
                freq4n.at(k0) = exp(LR_b(k*typenum+4));
                freq5n.at(k0) = exp(LR_b(k*typenum+5));
                freq6n.at(k0) = exp(LR_b(k*typenum+6));
                freq7n.at(k0) = exp(LR_b(k*typenum+7));
                freq8n.at(k0) = exp(LR_b(k*typenum+8));
                freq9n.at(k0) = exp(LR_b(k*typenum+9));
                freq10n.at(k0) = exp(LR_b(k*typenum+10));
                freq11n.at(k0) = exp(LR_b(k*typenum+11));
                freq12n.at(k0) = exp(LR_b(k*typenum+12));
                freq13n.at(k0) = exp(LR_b(k*typenum+13));
                freq14n.at(k0) = exp(LR_b(k*typenum+14));
                double sum = freq0n.at(k0)+freq1n.at(k0)+freq2n.at(k0)+freq3n.at(k0)+freq4n.at(k0)+freq5n.at(k0)+freq6n.at(k0)+freq7n.at(k0)+freq8n.at(k0)+freq9n.at(k0)+freq10n.at(k0)+freq11n.at(k0)+freq11n.at(k0)+freq12n.at(k0)+freq13n.at(k0)+freq14n.at(k0)+1;
                freq0n.at(k0) = freq0n.at(k0)/sum;
                freq1n.at(k0) = freq1n.at(k0)/sum;
                freq2n.at(k0) = freq2n.at(k0)/sum;
                freq3n.at(k0) = freq3n.at(k0)/sum;
                freq4n.at(k0) = freq4n.at(k0)/sum;
                freq5n.at(k0) = freq5n.at(k0)/sum;
                freq6n.at(k0) = freq6n.at(k0)/sum;
                freq7n.at(k0) = freq7n.at(k0)/sum;
                freq8n.at(k0) = freq8n.at(k0)/sum;
                freq9n.at(k0) = freq9n.at(k0)/sum;
                freq10n.at(k0) = freq10n.at(k0)/sum;
                freq11n.at(k0) = freq11n.at(k0)/sum;
                freq12n.at(k0) = freq12n.at(k0)/sum;
                freq13n.at(k0) = freq13n.at(k0)/sum;
                freq14n.at(k0) = freq14n.at(k0)/sum;
                freq15n.at(k0) = (double)1/sum;
            }
            Results[i].LR_freq_n.push_back(freq15n);
            Results[i].LR_freq_n.push_back(freq14n);
            Results[i].LR_freq_n.push_back(freq13n);
            Results[i].LR_freq_n.push_back(freq12n);
            Results[i].LR_freq_n.push_back(freq11n);
            Results[i].LR_freq_n.push_back(freq10n);
            Results[i].LR_freq_n.push_back(freq9n);
            Results[i].LR_freq_n.push_back(freq8n);
            Results[i].LR_freq_n.push_back(freq7n);
            Results[i].LR_freq_n.push_back(freq6n);
            Results[i].LR_freq_n.push_back(freq5n);
            Results[i].LR_freq_n.push_back(freq4n);
            Results[i].LR_freq_n.push_back(freq3n);
            Results[i].LR_freq_n.push_back(freq2n);
            Results[i].LR_freq_n.push_back(freq1n);
            Results[i].LR_freq_n.push_back(freq0n);
            
            VectorXd MLR_b = X*MLR_coeff;
            for (int k=0; k<numppos; k++){
                freq0p.at(numppos-k-1) = exp(MLR_b(k*typenum));
                freq1p.at(numppos-k-1) = exp(MLR_b(k*typenum+1));
                freq2p.at(numppos-k-1) = exp(MLR_b(k*typenum+2));
                freq3p.at(numppos-k-1) = exp(MLR_b(k*typenum+3));
                freq4p.at(numppos-k-1) = exp(MLR_b(k*typenum+4));
                freq5p.at(numppos-k-1) = exp(MLR_b(k*typenum+5));
                freq6p.at(numppos-k-1) = exp(MLR_b(k*typenum+6));
                freq7p.at(numppos-k-1) = exp(MLR_b(k*typenum+7));
                freq8p.at(numppos-k-1) = exp(MLR_b(k*typenum+8));
                freq9p.at(numppos-k-1) = exp(MLR_b(k*typenum+9));
                freq10p.at(numppos-k-1) = exp(MLR_b(k*typenum+10));
                freq11p.at(numppos-k-1) = exp(MLR_b(k*typenum+11));
                freq12p.at(numppos-k-1) = exp(MLR_b(k*typenum+12));
                freq13p.at(numppos-k-1) = exp(MLR_b(k*typenum+13));
                freq14p.at(numppos-k-1) = exp(MLR_b(k*typenum+14));
                double sum = freq0p.at(numppos-k-1)+freq1p.at(numppos-k-1)+freq2p.at(numppos-k-1)+freq3p.at(numppos-k-1)+freq4p.at(numppos-k-1)+freq5p.at(numppos-k-1)+freq6p.at(numppos-k-1)+freq7p.at(numppos-k-1)+freq8p.at(numppos-k-1)+freq9p.at(numppos-k-1)+freq10p.at(numppos-k-1)+freq11p.at(numppos-k-1)+freq11p.at(numppos-k-1)+freq12p.at(numppos-k-1)+freq13p.at(numppos-k-1)+freq14p.at(numppos-k-1)+1;
                freq0p.at(numppos-k-1) = freq0p.at(numppos-k-1)/sum;
                freq1p.at(numppos-k-1) = freq1p.at(numppos-k-1)/sum;
                freq2p.at(numppos-k-1) = freq2p.at(numppos-k-1)/sum;
                freq3p.at(numppos-k-1) = freq3p.at(numppos-k-1)/sum;
                freq4p.at(numppos-k-1) = freq4p.at(numppos-k-1)/sum;
                freq5p.at(numppos-k-1) = freq5p.at(numppos-k-1)/sum;
                freq6p.at(numppos-k-1) = freq6p.at(numppos-k-1)/sum;
                freq7p.at(numppos-k-1) = freq7p.at(numppos-k-1)/sum;
                freq8p.at(numppos-k-1) = freq8p.at(numppos-k-1)/sum;
                freq9p.at(numppos-k-1) = freq9p.at(numppos-k-1)/sum;
                freq10p.at(numppos-k-1) = freq10p.at(numppos-k-1)/sum;
                freq11p.at(numppos-k-1) = freq11p.at(numppos-k-1)/sum;
                freq12p.at(numppos-k-1) = freq12p.at(numppos-k-1)/sum;
                freq13p.at(numppos-k-1) = freq13p.at(numppos-k-1)/sum;
                freq14p.at(numppos-k-1) = freq14p.at(numppos-k-1)/sum;
                freq15p.at(numppos-k-1) = (double)1/sum;
            }
            Results[i].MLR_freq_p.push_back(freq0p);
            Results[i].MLR_freq_p.push_back(freq1p);
            Results[i].MLR_freq_p.push_back(freq2p);
            Results[i].MLR_freq_p.push_back(freq3p);
            Results[i].MLR_freq_p.push_back(freq4p);
            Results[i].MLR_freq_p.push_back(freq5p);
            Results[i].MLR_freq_p.push_back(freq6p);
            Results[i].MLR_freq_p.push_back(freq7p);
            Results[i].MLR_freq_p.push_back(freq8p);
            Results[i].MLR_freq_p.push_back(freq9p);
            Results[i].MLR_freq_p.push_back(freq10p);
            Results[i].MLR_freq_p.push_back(freq11p);
            Results[i].MLR_freq_p.push_back(freq12p);
            Results[i].MLR_freq_p.push_back(freq13p);
            Results[i].MLR_freq_p.push_back(freq14p);
            Results[i].MLR_freq_p.push_back(freq15p);

            for (int k=numppos; k<numnpos+numnpos; k++){
                int k0 = k-numppos;
                freq0n.at(k0) = exp(MLR_b(k*typenum));
                freq1n.at(k0) = exp(MLR_b(k*typenum+1));
                freq2n.at(k0) = exp(MLR_b(k*typenum+2));
                freq3n.at(k0) = exp(MLR_b(k*typenum+3));
                freq4n.at(k0) = exp(MLR_b(k*typenum+4));
                freq5n.at(k0) = exp(MLR_b(k*typenum+5));
                freq6n.at(k0) = exp(MLR_b(k*typenum+6));
                freq7n.at(k0) = exp(MLR_b(k*typenum+7));
                freq8n.at(k0) = exp(MLR_b(k*typenum+8));
                freq9n.at(k0) = exp(MLR_b(k*typenum+9));
                freq10n.at(k0) = exp(MLR_b(k*typenum+10));
                freq11n.at(k0) = exp(MLR_b(k*typenum+11));
                freq12n.at(k0) = exp(MLR_b(k*typenum+12));
                freq13n.at(k0) = exp(MLR_b(k*typenum+13));
                freq14n.at(k0) = exp(MLR_b(k*typenum+14));
                double sum = freq0n.at(k0)+freq1n.at(k0)+freq2n.at(k0)+freq3n.at(k0)+freq4n.at(k0)+freq5n.at(k0)+freq6n.at(k0)+freq7n.at(k0)+freq8n.at(k0)+freq9n.at(k0)+freq10n.at(k0)+freq11n.at(k0)+freq11n.at(k0)+freq12n.at(k0)+freq13n.at(k0)+freq14n.at(k0)+1;
                freq0n.at(k0) = freq0n.at(k0)/sum;
                freq1n.at(k0) = freq1n.at(k0)/sum;
                freq2n.at(k0) = freq2n.at(k0)/sum;
                freq3n.at(k0) = freq3n.at(k0)/sum;
                freq4n.at(k0) = freq4n.at(k0)/sum;
                freq5n.at(k0) = freq5n.at(k0)/sum;
                freq6n.at(k0) = freq6n.at(k0)/sum;
                freq7n.at(k0) = freq7n.at(k0)/sum;
                freq8n.at(k0) = freq8n.at(k0)/sum;
                freq9n.at(k0) = freq9n.at(k0)/sum;
                freq10n.at(k0) = freq10n.at(k0)/sum;
                freq11n.at(k0) = freq11n.at(k0)/sum;
                freq12n.at(k0) = freq12n.at(k0)/sum;
                freq13n.at(k0) = freq13n.at(k0)/sum;
                freq14n.at(k0) = freq14n.at(k0)/sum;
                freq15n.at(k0) = (double)1/sum;
            }
            Results[i].MLR_freq_n.push_back(freq15n);
            Results[i].MLR_freq_n.push_back(freq14n);
            Results[i].MLR_freq_n.push_back(freq13n);
            Results[i].MLR_freq_n.push_back(freq12n);
            Results[i].MLR_freq_n.push_back(freq11n);
            Results[i].MLR_freq_n.push_back(freq10n);
            Results[i].MLR_freq_n.push_back(freq9n);
            Results[i].MLR_freq_n.push_back(freq8n);
            Results[i].MLR_freq_n.push_back(freq7n);
            Results[i].MLR_freq_n.push_back(freq6n);
            Results[i].MLR_freq_n.push_back(freq5n);
            Results[i].MLR_freq_n.push_back(freq4n);
            Results[i].MLR_freq_n.push_back(freq3n);
            Results[i].MLR_freq_n.push_back(freq2n);
            Results[i].MLR_freq_n.push_back(freq1n);
            Results[i].MLR_freq_n.push_back(freq0n);
        }
    }
}

void freqplt(vector<res > &Results, size_t ** Table, int model, int numppos, int numnpos, int order, double like, string &figname){
    plt::figure_size(1200, 780);
    plt::suptitle(figname+", order is "+to_string(order)+", total log-likelihood is "+to_string(like));
    for (int k=0;k<Results.size();k++){
        int k0;
        if (model < 2){
            k0 = floor(k/2);
        }else{
            k0 = k;
        }
        int m;
        if (Results[k].LR_freq_p.size()>0){
            m = Results[k].LR_freq_p.size();
        }else{
            m = Results[k].LR_freq_n.size();
        }
        for (int j=0;j<m;j++){
            plt::subplot(4, 4, j+1+4*k0);
            //cout<<"Figure "<<j+1+4*k0<<"\n";
            if (Results[k].dir==0){
                int n = numppos+numnpos;
                std::vector<double> x(n), freq(n), LR_freq(n), MLR_freq(n);
                for (int i=0;i<numnpos;i++){
                    x.at(i) = -numnpos+i;
                    if (model%2==1){
                        freq.at(i) = (double)Table[numppos+i][4+j+4*k0]/(double)(Table[numppos+i][4+4*k0]+Table[numppos+i][5+4*k0]+Table[numppos+i][6+4*k0]+Table[numppos+i][7+4*k0]);
                    }else{
                        freq.at(i) = (double)Table[numppos+i][4+j+4*k0]/(double)(Table[numppos+i][4]+Table[numppos+i][5]+Table[numppos+i][6]+Table[numppos+i][7]+Table[numppos+i][8]+Table[numppos+i][9]+Table[numppos+i][10]+Table[numppos+i][11]+Table[numppos+i][12]+Table[numppos+i][13]+Table[numppos+i][14]+Table[numppos+i][15]+Table[numppos+i][16]+Table[numppos+i][17]+Table[numppos+i][18]+Table[numppos+i][19]);
                    }
                    LR_freq.at(i) = Results[k].LR_freq_n[j][i];
                    MLR_freq.at(i) = Results[k].MLR_freq_n[j][i];
                }
                for (int i=numnpos;i<numnpos+numppos;i++){
                    x.at(i) = i-numnpos+1;
                    if (model%2==1){
                        freq.at(i) = (double)Table[numnpos+numppos-1-i][4+j+4*k0]/(double)(Table[numnpos+numppos-1-i][4+4*k0]+Table[numnpos+numppos-1-i][5+4*k0]+Table[numnpos+numppos-1-i][6+4*k0]+Table[numnpos+numppos-1-i][7+4*k0]);
                    }else{
                        freq.at(i) = (double)Table[numnpos+numppos-1-i][4+j+4*k0]/(double)(Table[numnpos+numppos-1-i][4]+Table[numnpos+numppos-1-i][5]+Table[numnpos+numppos-1-i][6]+Table[numnpos+numppos-1-i][7]+Table[numnpos+numppos-1-i][8]+Table[numnpos+numppos-1-i][9]+Table[numnpos+numppos-1-i][10]+Table[numnpos+numppos-1-i][11]+Table[numnpos+numppos-1-i][12]+Table[numnpos+numppos-1-i][13]+Table[numnpos+numppos-1-i][14]+Table[numnpos+numppos-1-i][15]+Table[numnpos+numppos-1-i][16]+Table[numnpos+numppos-1-i][17]+Table[numnpos+numppos-1-i][18]+Table[numnpos+numppos-1-i][19]);
                    }
                    LR_freq.at(i) = Results[k].LR_freq_p[j][i-numnpos];
                    MLR_freq.at(i) = Results[k].MLR_freq_p[j][i-numnpos];
                }
                plt::plot(x,freq,"xr");
                plt::plot(x,LR_freq,"-ob");
                plt::plot(x,MLR_freq,"-om");
            }else if(Results[k].dir<0){
                //cout<<"Negative "<<k0<<"\n";
                int n = numnpos;
                std::vector<double> x(n), freq(n), LR_freq(n), MLR_freq(n);
                for (int i=0;i<numnpos;i++){
                    x.at(i) = -numnpos+i;
                    if (model%2==1){
                        freq.at(i) = (double)Table[numppos+i][4+j+4*k0]/(double)(Table[numppos+i][4+4*k0]+Table[numppos+i][5+4*k0]+Table[numppos+i][6+4*k0]+Table[numppos+i][7+4*k0]);
                    }else{
                        freq.at(i) = (double)Table[numppos+i][4+j+4*k0]/(double)(Table[numppos+i][4]+Table[numppos+i][5]+Table[numppos+i][6]+Table[numppos+i][7]+Table[numppos+i][8]+Table[numppos+i][9]+Table[numppos+i][10]+Table[numppos+i][11]+Table[numppos+i][12]+Table[numppos+i][13]+Table[numppos+i][14]+Table[numppos+i][15]+Table[numppos+i][16]+Table[numppos+i][17]+Table[numppos+i][18]+Table[numppos+i][19]);
                    }
                    LR_freq.at(i) = Results[k].LR_freq_n[j][i];
                    MLR_freq.at(i) = Results[k].MLR_freq_n[j][i];
                }
                plt::plot(x,freq,"xr");
                plt::plot(x,LR_freq,"-ob");
                plt::plot(x,MLR_freq,"-om");
            }else{
                //cout<<"Positive "<<k0<<"\n";
                int n = numnpos;
                std::vector<double> x(n), freq(n), LR_freq(n), MLR_freq(n);
                for (int i=0;i<numppos;i++){
                    x.at(i) = i+1;
                    if (model%2==1){
                        freq.at(i) = (double)Table[numppos-1-i][4+j+4*k0]/(double)(Table[numppos-1-i][4+4*k0]+Table[numppos-1-i][5+4*k0]+Table[numppos-1-i][6+4*k0]+Table[numppos-1-i][7+4*k0]);
                    }else{
                        freq.at(i) = (double)Table[numppos-1-i][4+j+4*k0]/(double)(Table[numppos-1-i][4]+Table[numppos-1-i][5]+Table[numppos-1-i][6]+Table[numppos-1-i][7]+Table[numppos-1-i][8]+Table[numppos-1-i][9]+Table[numppos-1-i][10]+Table[numppos-1-i][11]+Table[numppos-1-i][12]+Table[numppos-1-i][13]+Table[numppos-1-i][14]+Table[numppos-1-i][15]+Table[numppos-1-i][16]+Table[numppos-1-i][17]+Table[numppos-1-i][18]+Table[numppos-1-i][19]);
                    }
                    LR_freq.at(i) = Results[k].LR_freq_p[j][i];
                    MLR_freq.at(i) = Results[k].MLR_freq_p[j][i];
                }
                plt::plot(x,freq,"xr");
                plt::plot(x,LR_freq,"-ob");
                plt::plot(x,MLR_freq,"-om");
            }
        }
    }
    plt::savefig("regressionfig.png");
    plt::show();
    
}

void output(const char* outputname, const char* filename, string figname, vector<res > &Results, int model, int order, double like){
      string ss = "ACGT";
      ofstream myfile;
      myfile.open(outputname);
      string s(filename);
      myfile << "Output of the regression of "+s+"\n";
      myfile <<"Method: "<< figname +"\n";
      myfile <<"Regression order is "<< to_string(order) +"\n";
      myfile <<"Total likelihood is "<< to_string(like) +"\n";
      myfile <<"\n";
    
    if (model<=1){
        for (int i=0; i<Results.size(); i++){
            if (Results[i].dir>0){
                myfile << "Direction: Forward Reference nucleotide "+Results[i].refnuc+", the regression results obey forward order.\n";
                myfile << "Logit regression results\tMultinomial logit regression results\n";
                for (int j=0; j<Results[i].LR_coeff.size(); j++){
                    myfile << to_string(Results[i].LR_coeff(j))<<"\t"<<to_string(Results[i].MLR_coeff(j))<<"\n";
                }
            }else{
                myfile << "Direction: Backward Reference nucleotide "+Results[i].refnuc+", the regression results obey backward order.\n";
                myfile << "Logit regression results\tMultinomial logit regression results\n";
                for (int j=0; j<Results[i].LR_coeff.size(); j++){
                    myfile << to_string(Results[i].LR_coeff(j))<<"\t"<<to_string(Results[i].MLR_coeff(j))<<"\n";
                }
            }
        }
    }else{
        for (int i=0; i<Results.size(); i++){
            myfile << "Direction: Forward Reference nucleotide "+Results[i].refnuc+"\t"+"Backward Reference nucleotide "+ss[3-Results[i].refnucid]+", the regression results obey forward order.\n";
            myfile << "Logit regression results\tMultinomial logit regression results\n";
            for (int j=0; j<Results[i].LR_coeff.size(); j++){
                myfile << to_string(Results[i].LR_coeff(j))<<"\t"<<to_string(Results[i].MLR_coeff(j))<<"\n";
            }
        }
    }
      
      myfile.close();
}

int main_regression(int argc,char**argv){
    int numppos = 15;
    int numnpos = 15;
    int numcolumn = 16 + 4;
    int numpos = numppos + numnpos;
    int order = 6;
    
    size_t ** Table = (size_t **) malloc(numpos*(sizeof(int *))); /*I allocate memory here.  If this function is called many times it may be better to move the memmory allocation out of this function*/
    for (int i=0; i<numpos; i++){
        Table[i]=(size_t *) malloc(numcolumn*(sizeof(size_t)));
    }
    string* ColumnName = (string *) malloc(numcolumn*(sizeof(string)));
    const char* filename = "data_ancient_human.txt";
    readdata(filename, ColumnName,Table);
    
    int model = 1;
    cout<< Table[0][0] <<"\n";
    cout<< numppos <<"\n";
    
    vector<res > Results;
    double like = 0;
    string figname;
    modelcalculation(Table, numppos, numnpos, model, order, Results, figname);
    for (int i = 0;i<Results.size(); i++){
        like = like + Results[i].MLR_like;
    }
    cout<<"Likelihood is "<<like<<"\n";
    if (model%2 == 0){
        Freqcalculator4uncon(Results, numppos, numnpos);
    }else{
        Freqcalculator4cond(Results, numppos, numnpos);
    }
    //freqplt(Results, Table, model, numppos, numnpos, order, like, figname);
    output("test.out", filename, figname, Results, model, order, like);
    
    for (int i = 0; i < numpos; i++)
    free(Table[i]);
    free(Table);
    
    return 0;
}

//int main(){
//    main_regression();
//}
