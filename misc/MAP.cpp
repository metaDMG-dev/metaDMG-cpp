#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <iostream>
#include <math.h>
#include <cassert>

#include "../bfgs.h"
#include "M3Load.h"
#include "Likelihood.h"

#define NUMROWS 30
#define KCOL 21
#define NCOL 22
#define XCOL 24
#define C_freq 19
#define G_freq 20

typedef struct{
    const char* M3_matrix_print;
    const char* OutputStat;
}argStruct;

int HelpPage(FILE *fp){
  fprintf(fp,"Generate artificial reference genome uniformly sampled from the four nucleotides\n");
  fprintf(fp,"Usage\n./MAP -m <Suffix.mismatches.txt.gz> -o <NAME.txt.gz>\n");
  fprintf(fp,"\nExample\n\n");
  fprintf(fp,"\nOptions: \n");
  fprintf(fp,"-h   | --help: \t\t\t Print help page.\n");
  fprintf(fp,"-v   | --version: \t\t Print help page.\n\n");
  fprintf(fp,"-m  | --mismatchmatrix: \t Mismatch matrix from metadamage print or print_ugly, not bdamage file\n");
  fprintf(fp,"-o   | --outstat: \t\t output statistics file in .gz format\n");
  exit(1);
  return 0;
}
 
argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->M3_matrix_print = NULL;
  mypars->OutputStat = NULL;
  ++argv;
  while(*argv){
    //fprintf(stderr,"ARGV %s\n",*argv);
    if(strcasecmp("-m",*argv)==0 || strcasecmp("--mismatchmatrix",*argv)==0){
      mypars->M3_matrix_print = strdup(*(++argv));
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--outstat",*argv)==0){
      mypars->OutputStat = strdup(*(++argv));
    }
    else{
      fprintf(stderr,"unrecognized input option %s, see help page\n\n",*(argv));
      exit(0);
    }
    ++argv;
  }
  return mypars;
}

void MAP(double M3[MAX_ROWS][MAX_COLS],double* MAP){
    double lowbound[] = {0.00000001,0.00000001,0.00000001,2};
    double upbound[] = {1-0.00000001,1-0.00000001,1-0.00000001,100000};
    int nbd[] = {2,2,2,1};
    int noisy = 0;
    int numpars = 4;//sizeof(MAP) / sizeof(double); // determine length of MAP array;
    double* DMGparam = (double*) malloc(numpars*sizeof(double)); // A q c phi

    DMGparam[0] = 0.1;
    DMGparam[1] = 0.1;
    DMGparam[2] = 0.01; 
    DMGparam[3] = 1000; 

    fprintf(stderr,"I AM RUNNING findmax_bfgs NOW\n");
    findmax_bfgs(numpars, DMGparam, nullptr,compute_log_likelihood,nullptr,lowbound,upbound,nbd,noisy);
    fprintf(stderr,"OPTIM IS DONE NOW\n");

    MAP[0] = DMGparam[0];
    MAP[1] = DMGparam[1];
    MAP[2] = DMGparam[2];
    MAP[3] = DMGparam[3];
    MAP[4] = compute_log_likelihood(DMGparam,nullptr);

    free(DMGparam);
}

void MAP_stat_file(double M3[MAX_ROWS][MAX_COLS],double* LlhRes,gzFile file){

    // MAP damage, damage_std Map_damage_significance
    int N_pos = M3[0][NCOL];

    // mu = A = LlhRes[0]
    double std = std::sqrt((LlhRes[0]*(1-LlhRes[0])*(LlhRes[3]+N_pos))/((LlhRes[3]+1)*N_pos));
    double significance = LlhRes[0]/std;

    // MAP_Damage   MAP_signifcance     MAP_damage_std
    gzprintf(file,"%f \t %f \t %f \t",LlhRes[0],std,significance);

    /*
    add the significance etc for the A, q, c when done
    */

    double q = LlhRes[1];
    double A = LlhRes[0];
    double c = LlhRes[2];
    double phi = LlhRes[3];
    double llh = LlhRes[4];
    // then the optimized values MAP_A MAP_q MAP_phi MAP_c 
    gzprintf(file,"%f \t %f \t %f \t %f \t %f \t",A,q,phi,c,llh);
    
    double alpha;
    double beta;
    double Dx_var_numerator;
    double Dx_var_deumerator;
    double Dx_std;
    double Dx_std_norm;
    for(int i = 0; i < NUMROWS;i++){
        double Dx = A * pow((1 - q), fabs(M3[i][XCOL]) - 1) + c;
        gzprintf(file,"%f \t",Dx);
    }

    for(int i = 0; i < NUMROWS;i++){
        double Dx = A * pow((1 - q), fabs(M3[i][XCOL]) - 1) + c;
        alpha = Dx * phi;
        beta = (1 - Dx) * phi;
        Dx_var_numerator = M3[i][NCOL]*alpha*beta*(alpha+beta+M3[i][NCOL]);
        Dx_var_deumerator = pow((alpha+beta),2)*(alpha+beta+1);
        Dx_std = std::sqrt((Dx_var_numerator/Dx_var_deumerator));
        Dx_std_norm = Dx_std / M3[i][NCOL];
        gzprintf(file,"%f",Dx_std_norm);
        
        if (i < NUMROWS-1){
            //ensures the elements are tab seperated
            gzprintf(file,"\t");
        }
        else{
            // except for the last making new line
            gzprintf(file,"\n");
        }
    }

    /*mu = A
    std = np.sqrt(A * (1 - A) * (phi + N) / ((phi + 1) * N))
    significance = mu / std
    return mu, std, significance
    
    A:  0.21290688706956776  q:  0.35714159577362814  c:  0.009605540202966379  phi:  9369.817278739834
N MAX IS  28921
 std is  0.00486586260361796
 Sig is  43.75521966264792
 
 */
}


void M3_stat_file(double M3[MAX_ROWS][MAX_COLS],gzFile file){
    
    int N_max = 0;
    
    // N_x=1_forward  N_x=1_reverse
    gzprintf(file,"%d \t %d \t",(int) M3[0][NCOL],(int)M3[NUMROWS-1][NCOL]);

    int N_sum = 0;
    int N_sum_fwd = 0;
    int N_sum_rev = 0;
    int K_sum_fwd = 0;
    int K_sum_rev = 0;
    for(int i = 0; i < NUMROWS;i++){
        N_sum += M3[i][NCOL];
        if (i < NUMROWS/2){
            N_sum_fwd += M3[i][NCOL];
            K_sum_fwd += M3[i][KCOL];
        }
        else if (i > NUMROWS/2){
            N_sum_rev += M3[i][NCOL];
            K_sum_rev += M3[i][KCOL];
        }
    }

    gzprintf(file,"%d \t %d \t %d \t %d \t %d\t",N_sum,N_sum_fwd,N_sum_rev,(int)M3[0][26],(int)M3[0][27]);
    gzprintf(file,"%d \t %d \t %d \t",(int)M3[0][25],K_sum_fwd,K_sum_rev);
    /*
            M3[i][25] = k_sum_total;
        M3[i][26] = Min_N_in_group;
        M3[i][27] = Max_N_in_group;
        */
    // 5'
    for(int i = 0; i < NUMROWS/2;i++){
        //fprintf(stderr,"K col %d \n",i);
        gzprintf(file,"%d \t",(int) M3[i][KCOL]);
    }

    for(int i = 0; i < NUMROWS/2;i++){
        //fprintf(stderr,"N col %d \n",i);
        gzprintf(file,"%d \t",(int) M3[i][NCOL]);
    }

    for(int i = 0; i < NUMROWS/2;i++){
        //fprintf(stderr,"F col %d \n",i);
        gzprintf(file,"%f \t",M3[i][C_freq]);
    }

    // 3'
    for(int i = NUMROWS; i > NUMROWS/2;i--){
        //fprintf(stderr,"K col %d \n",i);
        gzprintf(file,"%d \t",(int) M3[i-1][KCOL]);
    }

     for(int i = NUMROWS; i > NUMROWS/2;i--){
        //fprintf(stderr,"N col %d \n",i);
        gzprintf(file,"%d \t",(int) M3[i-1][NCOL]);
    }
    int last_i = 0;
    for(int i = NUMROWS; i > NUMROWS/2+1;i--){
        //fprintf(stderr,"f col %d \n",i);
        gzprintf(file,"%f \t",M3[i-1][G_freq]);
        last_i = i;
    }
    gzprintf(file,"%f \t",M3[last_i-1][G_freq]);
}

const char** getColumnNames(int* colnumber) {
    const int maxColumnNames = 200;
    const char** columnNames = (const char**)malloc(200 * sizeof(const char*));
    if (columnNames == NULL) {
        fprintf(stderr,"Failed to allocate memory for column names");
        return NULL;
    }

    columnNames[0] = "sample";
    columnNames[1] = "tax_id";

    int N_sum = 0;
    int N_sum_fwd = 0;
    int N_sum_rev = 0;
    int K_sum_fwd = 0;
    int K_sum_rev = 0;

    for (int i = 0; i < NUMROWS; i++) {
        N_sum += M3[i][NCOL];
        if (i < NUMROWS / 2) {
            N_sum_fwd += M3[i][NCOL];
            K_sum_fwd += M3[i][KCOL];
        }
        else if (i > NUMROWS / 2) {
            N_sum_rev += M3[i][NCOL];
            K_sum_rev += M3[i][KCOL];
        }
    }

    columnNames[2] = "N_x=1_forward";
    columnNames[3] = "N_x=1_reverse";
    columnNames[4] = "N_sum_total";
    columnNames[5] = "N_sum_forward";
    columnNames[6] = "N_sum_reverse";
    columnNames[7] = "N_min";
    columnNames[8] = "N_max";

    columnNames[9] = "k_sum";
    columnNames[10] = "k_sum_forward";
    columnNames[11] = "k_sum_reverse";

    int curr_it = 11;

    for (int i = 0; i < NUMROWS / 2; i++) {
        char k_buff_fwd[10];
        snprintf(k_buff_fwd, sizeof(k_buff_fwd), "k+%d", i + 1);
        columnNames[curr_it+i+1] = strdup(k_buff_fwd);
    }

    curr_it += NUMROWS / 2;
    for (int i = 0; i < NUMROWS / 2; i++) {
        char N_buff_fwd[10];
        snprintf(N_buff_fwd, sizeof(N_buff_fwd), "N+%d", i + 1);
        columnNames[curr_it+i+1] = strdup(N_buff_fwd);
    }

    curr_it += NUMROWS / 2;
    for (int i = 0; i < NUMROWS / 2; i++) {
        char f_buff_fwd[10];
        snprintf(f_buff_fwd, sizeof(f_buff_fwd), "f+%d", i + 1);
        columnNames[curr_it+i+1] = strdup(f_buff_fwd);
    }

    curr_it += NUMROWS / 2;
    for (int i = 0; i < NUMROWS / 2; i++) {
        char k_buff_rev[10];
        snprintf(k_buff_rev, sizeof(k_buff_rev), "k-%d", i + 1);
        columnNames[curr_it+i+1] = strdup(k_buff_rev);
    }

    curr_it += NUMROWS / 2;
    for (int i = 0; i < NUMROWS / 2; i++) {
        char N_buff_rev[10];
        snprintf(N_buff_rev, sizeof(N_buff_rev), "N-%d", i + 1);
        columnNames[curr_it+i+1] = strdup(N_buff_rev);
    }
    
    curr_it += NUMROWS / 2;
    for (int i = 0; i < NUMROWS / 2; i++) {
        char f_buff_rev[10];
        snprintf(f_buff_rev, sizeof(f_buff_rev), "f-%d", i + 1);
        columnNames[curr_it+i+1] = strdup(f_buff_rev);
    }
    
    curr_it += NUMROWS / 2;
    columnNames[curr_it + 1] = "MAP";
    columnNames[curr_it + 2] = "MAP_damage_std";
    columnNames[curr_it + 3] = "MAP_damage_significance";
    columnNames[curr_it + 4] = "A";
    columnNames[curr_it + 5] = "q";
    columnNames[curr_it + 6] = "phi";
    columnNames[curr_it + 7] = "c";
    columnNames[curr_it + 8] = "llh";

    curr_it += 8;
    for (int i = 0; i < NUMROWS / 2; i++) {
        //std::cout << "values curr_it " << curr_it + i << std::endl;
        char Dx_buff_fwd[10];
        snprintf(Dx_buff_fwd, sizeof(Dx_buff_fwd), "Dx+%d", i + 1);
        columnNames[curr_it+i+1] = strdup(Dx_buff_fwd);
    }

    curr_it += NUMROWS / 2;
    for (int i = 0; i < NUMROWS/2; i++) {
        char Dx_buff_rev[10];
        snprintf(Dx_buff_rev,10, "Dx-%d", i + 1);
        columnNames[curr_it+i+1] = strdup(Dx_buff_rev);
    }
    
    curr_it += NUMROWS / 2;
    for (int i = 0; i < NUMROWS/2; i++) {
        char Dx_std_buff_fwd[10];
        snprintf(Dx_std_buff_fwd,10, "Dx_std+%d", i + 1);
        columnNames[curr_it+i+1] = strdup(Dx_std_buff_fwd);
    }

    curr_it += NUMROWS / 2;
    for (int i = 0; i < NUMROWS/2; i++) {
        char Dx_std_buff_rev[10];
        snprintf(Dx_std_buff_rev,10, "Dx_std-%d", i + 1);
        columnNames[curr_it+i+1] = strdup(Dx_std_buff_rev);
    }
    curr_it += NUMROWS / 2;

    (*colnumber) = curr_it+1;
    return columnNames;
}

void M3Print_to_OutStat(int argc,char **argv){
    argStruct *mypars = NULL;
    if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"-v")==0||
                            strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
        HelpPage(stderr);
    }
    else{
        mypars = getpars(argc,argv);

        int num_rows, num_cols;
        const char* M3file = mypars-> M3_matrix_print; //"MycoBactBamSEOutSortMDSortN.mismatches.txt.gz";
        const char* filename = mypars-> OutputStat; //"outputtest.txt.gz";

        read_count_matrix(M3file, M3,tax_id,dir,&num_rows, &num_cols);
        Alter_count_matrix(M3,tax_id,dir,num_rows,num_cols);

        int numpars = 5;
        double* LlhRes = (double*) malloc(numpars*sizeof(double));    
        MAP(M3,LlhRes);
        //fprintf(stderr,"A: %f \t q: %f \t c: %f \t phi: %f \t llh: %f \n",LlhRes[0],LlhRes[1],LlhRes[2],LlhRes[3],LlhRes[4]);

        int colnumber = 0;
        const char** columnNames = getColumnNames(&colnumber);

        gzFile gz = Z_NULL;
        gz = gzopen(filename,"w");
        assert(gz!=Z_NULL);

        if (columnNames != NULL) {
            for (int i = 0; i < colnumber; i++) {
                //fprintf(stderr,"%s \t",columnNames[i]);
                gzprintf(gz,"%s \t",columnNames[i]);
            }
            gzprintf(gz,"\n");    
            // Free the allocated memory
            free(columnNames);
        }

        char* Id_copy = strdup(M3file);
        char* sample_id = strtok(Id_copy,".");
        free(Id_copy);


        gzprintf(gz,"%s \t %s \t","tmp",tax_id[0]);
        M3_stat_file(M3,gz);
        MAP_stat_file(M3,LlhRes,gz);
        gzclose(gz);
    }
}




#ifdef __WITH_MAIN__
int main(int argc,char **argv){
    M3Print_to_OutStat(argc,argv); 
    return 0;
}
#endif

//g++ MAP.cpp M3Load.o Likelihood.o bfgs.o -std=c++11 -lm -lz -D __WITH_MAIN__ -o MAP 
