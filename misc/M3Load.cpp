#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <iostream>
#include "M3Load.h"

double M3[MAX_ROWS][MAX_COLS];
char tax_id[MAX_ROWS][MAX_COLS];
char dir[MAX_ROWS][MAX_COLS];

// Define the columns present in our M3 File after alterations
#define AA 1
#define AC 2
#define AG 3
#define AT 4

#define CA 5
#define CC 6
#define CG 7
#define CT 8

#define GA 9
#define GC 10
#define GG 11
#define GT 12

#define TA 13
#define TC 14
#define TG 15
#define TT 16

#define Ctotal 17
#define Gtotal 18
#define C_freq 19
#define G_freq 20

#define KCOL 21
#define NCOL 22
#define XCOL 24

void read_count_matrix(const char* filename, double M3[MAX_ROWS][MAX_COLS],char tax_id[MAX_ROWS][MAX_COLS],char dir[MAX_ROWS][MAX_COLS],int* num_rows, int* num_cols) {
    // Opens the mismatchmatrix file in .txt.gz format
    gzFile file = gzopen(filename, "rb");
    if (!file) {
        printf("Failed to open file: %s\n", filename);
        exit(1);
    }

    // Read header line and determine number of columns
    // ['tax_id', 'direction', 'position', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    char line[MAX_LINE_LENGTH];
    if (!gzgets(file, line, MAX_LINE_LENGTH)) {
        printf("Failed to read header line\n");
        exit(1);
    }
    
    char* token = strtok(line, "\t"); // Skip first column #taxidStr
    token = strtok(NULL, "\t"); // Skip second column, #direction

    // Read data lines and store counts in matrix to generate column dimension of M3
    *num_cols = 0;
    *num_rows = 0;
    while ((token = strtok(NULL, "\t")) != NULL) {
        (*num_cols)++;
    }

    while (gzgets(file, line, MAX_LINE_LENGTH)) {
        token = strtok(line, "\t"); // Skip first column #taxidStr
        strcpy(tax_id[*num_rows], token); 
        //fprintf(stderr,"THE TAXID %s FOR ROW %d IS \n",tax_id[*num_rows],*num_rows);
        token = strtok(NULL, "\t"); // Skip second column, direction
        strcpy(dir[*num_rows], token); 

        int col = 0;
        while ((token = strtok(NULL, "\t")) != NULL && col < *num_cols) {
            //fprintf(stderr," while look token %s \n",token);
            M3[*num_rows][col++] = atof(token);
        }
        (*num_rows)++;

    }
    gzclose(file);
}
    
void Alter_count_matrix(double M3[MAX_ROWS][MAX_COLS],char dir[MAX_ROWS][MAX_COLS],int num_rows){
    // Alters / adds additional columns to the count M^3 matrix, needed for log-likelihood calculations
    
    double C_total, G_total = 0; 

    int k_sum_total = 0;
    //    double N_in_group[MAX_ROWS];
    int Max_N_in_group = 0;
    int Min_N_in_group = 0;

    for (int i = 0; i < num_rows; i++){
        M3[i][XCOL] = abs(M3[i][0])+1; //df["|x|"] = np.abs(df["position"])
        M3[i][0]=i+1; //change the index of position
        C_total = M3[i][CA] + M3[i][CC] + M3[i][CG] + M3[i][CT];
        G_total = M3[i][GA] + M3[i][GC] + M3[i][GG] + M3[i][GT];

        M3[i][Ctotal] = C_total;
        M3[i][Gtotal] = G_total;
        M3[i][C_freq] = M3[i][CT]/C_total; //CT / C i.e freq C>T
        M3[i][G_freq] = M3[i][GA]/G_total; //GA / C i.e freq G>A
        if (strcasecmp(dir[i], "5'") == 0){
            M3[i][KCOL] = M3[i][CT]; // df["k"] = df["CT"]
            M3[i][NCOL] = C_total; //df["N"] = df["C"]
            M3[i][23] = M3[i][CT]/C_total; //df["f"] = df["k"] / df["N"]
        }
        else if (strcasecmp(dir[i], "3'") == 0){
            M3[i][KCOL] = M3[i][GA]; // df["k"] = df["CT"]
            M3[i][NCOL] = G_total; //df["N"] = df["C"]
            M3[i][23] = M3[i][GA]/G_total; //df["f"] = df["k"] / df["N"]
        }
        k_sum_total += (int)M3[i][21];
        
        if ((int)M3[i][NCOL] > Max_N_in_group){
            Max_N_in_group = (int) M3[i][NCOL];
            Min_N_in_group = (int)M3[i][NCOL]; //once the minimum is initated to the highest it is easier to decrease again
        }

        if ((int)M3[i][NCOL] < Min_N_in_group){
            Min_N_in_group = (int)M3[i][NCOL];
        }
    }
    
    //add additional columns
    for (int i = 0; i < num_rows; i++){
        M3[i][25] = 
        k_sum_total;
        M3[i][26] = Min_N_in_group;
        M3[i][27] = Max_N_in_group;
        
        /*
        fprintf(stderr,"dir %s \t \n%.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.4f %.4f %.4f %.4f %.4f %.4f %d %d %d\n",dir[i],
        M3[i][0],M3[i][1],M3[i][2],M3[i][3],M3[i][4],M3[i][5],M3[i][6],M3[i][7],
        M3[i][8],M3[i][9],M3[i][10],M3[i][11],M3[i][12],M3[i][13],M3[i][14],M3[i][15],
        M3[i][16],M3[i][17],M3[i][18],M3[i][19],M3[i][20],M3[i][21],M3[i][22],M3[i][23],M3[i][24],(int) M3[i][25],(int) M3[i][26],(int) M3[i][27]);
        */
    }

    // NOTE WE NEED A COLUMN 28 WHICH CONTAINS THE FILE NAME!
}

#ifdef __WITH_MAIN__
int main() {
    //double M3[MAX_ROWS][MAX_COLS];
    //char tax_id[MAX_ROWS][MAX_COLS];
    //char dir[MAX_ROWS][MAX_COLS];
    int num_rows, num_cols;
    read_count_matrix("MycoBactBamSEOutSortMDSortN.mismatches.txt.gz", M3,tax_id,dir,&num_rows, &num_cols);
    printf("Read %d rows and %d columns\n", num_rows, num_cols);
    Alter_count_matrix(M3,dir,num_rows);

    return 0;
}
#endif

//g++ M3Load.cpp -std=c++11 -lm -lz -o M3
//g++ M3Load.cpp -c -std=c++11 -lm -lz
//g++ M3Load.cpp -std=c++11 -lm -lz -D __WITH_MAIN__ -o M3
