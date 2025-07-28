#ifndef M3LOAD_H
#define M3LOAD_H
#define MAX_ROWS 100
#define MAX_COLS 100
#define MAX_LINE_LENGTH 1000
extern double M3[MAX_ROWS][MAX_COLS];
extern char tax_id[MAX_ROWS][MAX_COLS];
extern char dir[MAX_ROWS][MAX_COLS];

void read_count_matrix(const char* filename, double M3[MAX_ROWS][MAX_COLS],char tax_id[MAX_ROWS][MAX_COLS],char dir[MAX_ROWS][MAX_COLS],int* num_rows, int* num_cols);

void Alter_count_matrix(double M3[MAX_ROWS][MAX_COLS],char dir[MAX_ROWS][MAX_COLS],int num_rows);

#endif
