#ifndef PCA_H
#define PCA_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#define MAX_DIM 1024









void center_data(double *data, int m, int dim);








void compute_covariance_matrix(double *data, double *cov_matrix, int m, int dim);








void eigen_decompose(double *cov_matrix, double *eigenvalues, double *eigenvectors, int dim);












void project_data(double *data, double *eigenvectors, double *projected_data,
                  int m, int dim, int k, int *selected_indices, int output_samples);








int* generate_sample_indices(int total_samples, int output_samples, int sampling_method);











void PCA_transform_advanced(double **input, double **output, int m, int n, int k,
                            int output_samples, int sampling_method);









void PCA_transform(double **input, double **output, int m, int n, int k);







double** create_output_matrix(int rows, int cols);






void free_matrix(double **matrix, int rows);

#endif