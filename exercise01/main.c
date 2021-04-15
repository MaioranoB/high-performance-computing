#include <stdlib.h>
#include <stdio.h>
#include <time.h>

double *create_rand_array(int n);
double *create_empty_array(int n);
double **create_rand_square_matrix(int n);
double time_to_multiply_ij(int n, double **matrix, double *array, double *result);
double time_to_multiply_ji(int n, double **matrix, double *array, double *result);
void clear_array(int n, double *array);

int main(void){
    srand(time(NULL));

    int nMAX = 16000;
    
    // b = A * x
    double **A = create_rand_square_matrix(nMAX);
    double *x  = create_rand_array(nMAX);
    double *b  = create_empty_array(nMAX);

    FILE *csv_file = fopen("c_times.csv", "w");
    fprintf(csv_file,"n,Time_ij,Time_ji\n");

    int n;
    int range = 16;
    for (int i = 0; i <= range; i++){
        n = i * 1000;
        
        double time_ij = time_to_multiply_ij(n,A,x,b);
        clear_array(n,b);
        double time_ji = time_to_multiply_ji(n,A,x,b);
        
        printf("n = %5d -> ij: %.6f / ji: %.6f \n",n,time_ij,time_ji);
        fprintf(csv_file,"%d,%.6f,%.6f\n",n, time_ij, time_ji);
    }
    fclose(csv_file);
    return 0;
}

double *create_rand_array(int n){
    double *array = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++){
        array[i] = (double)rand();
    }
    return array;
}
double *create_empty_array(int n){
    double *array = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++){
        array[i] = 0;
    }
    return array;
}
double **create_rand_square_matrix(int n){
    double **matrix = (double **)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++){
        matrix[i] = create_rand_array(n);
    }
    return matrix;
}

double time_to_multiply_ij(int n, double **matrix, double *array, double *result){
    double s;
    clock_t start = clock();
    for (int i = 0; i < n; i++){
        s = 0;
        for (int j = 0; j < n; j++){
            // result[i] += matrix[i][j] * array[j];
            s += matrix[i][j] * array[j];
        }
        result[i] = s;
    }
    clock_t end = clock();
    return (double)(end - start)/CLOCKS_PER_SEC;
}
double time_to_multiply_ji(int n, double **matrix, double *array, double *result){
    double t;
    clock_t start = clock();
    for (int j = 0; j < n; j++){
        t = array[j];
        for (int i = 0; i < n; i++){
            // result[i] += matrix[i][j] * array[j];
            result[i] += matrix[i][j] * t;
        }
    }
    clock_t end = clock();
    return (double)(end - start)/CLOCKS_PER_SEC;
}

void clear_array(int n, double *array){
    for (int i = 0; i < n; i++){
        array[i] = 0;
    }
}
