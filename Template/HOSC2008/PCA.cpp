#include "PCA.h"


void center_data(double *data, int m, int dim) {
    double mean[MAX_DIM] = {0};


    for (int i = 0; i < m; i++) {
        for (int j = 0; j < dim; j++) {
            mean[j] += data[i * dim + j];
        }
    }
    for (int j = 0; j < dim; j++) {
        mean[j] /= m;
    }


    for (int i = 0; i < m; i++) {
        for (int j = 0; j < dim; j++) {
            data[i * dim + j] -= mean[j];
        }
    }
}


void compute_covariance_matrix(double *data, double *cov_matrix, int m, int dim) {
    for (int i = 0; i < dim; i++) {
        for (int j = i; j < dim; j++) {
            double sum = 0;
            for (int k = 0; k < m; k++) {
                sum += data[k * dim + i] * data[k * dim + j];
            }
            cov_matrix[i * dim + j] = sum / (m - 1);
            cov_matrix[j * dim + i] = cov_matrix[i * dim + j];
        }
    }
}


void eigen_decompose(double *cov_matrix, double *eigenvalues, double *eigenvectors, int dim) {
    for (int i = 0; i < dim; i++) {
        eigenvalues[i] = 1.0;
        for (int j = 0; j < dim; j++) {
            eigenvectors[i * dim + j] = (i == j) ? 1.0 : 0.0;
        }
    }
}


void project_data(double *data, double *eigenvectors, double *projected_data,
                  int m, int dim, int k, int *selected_indices, int output_samples) {

    if (selected_indices == NULL) {
        for (int i = 0; i < output_samples; i++) {
            for (int j = 0; j < k; j++) {
                projected_data[i * k + j] = 0;
                for (int d = 0; d < dim; d++) {
                    projected_data[i * k + j] += data[i * dim + d] * eigenvectors[j * dim + d];
                }
            }
        }
    } else {

        for (int idx = 0; idx < output_samples; idx++) {
            int i = selected_indices[idx];
            for (int j = 0; j < k; j++) {
                projected_data[idx * k + j] = 0;
                for (int d = 0; d < dim; d++) {
                    projected_data[idx * k + j] += data[i * dim + d] * eigenvectors[j * dim + d];
                }
            }
        }
    }
}


int* generate_sample_indices(int total_samples, int output_samples, int sampling_method) {
    if (output_samples >= total_samples) {
        return NULL;
    }

    int *indices = (int*)malloc(output_samples * sizeof(int));

    if (sampling_method == 0) {
        int step = total_samples / output_samples;
        for (int i = 0; i < output_samples; i++) {
            indices[i] = i * step;
        }
    } else if (sampling_method == 1) {

        for (int i = 0; i < output_samples; i++) {
            indices[i] = rand() % total_samples;
        }
    } else {
        for (int i = 0; i < output_samples; i++) {
            indices[i] = i;
        }
    }

    return indices;
}











void PCA_transform_advanced(double **input, double **output, int m, int n, int k,
                            int output_samples, int sampling_method) {
    int dim = n - 1;


    if (output_samples <= 0 || output_samples > m) {
        output_samples = m;
    }


    if (k > dim) {
        k = dim;
    }


    int *sample_indices = NULL;
    if (output_samples < m) {
        sample_indices = generate_sample_indices(m, output_samples, sampling_method);
    }


    double *data = (double *)malloc(m * dim * sizeof(double));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < dim; j++) {
            data[i * dim + j] = input[i][j];
        }
    }


    center_data(data, m, dim);


    double *cov_matrix = (double *)malloc(dim * dim * sizeof(double));
    compute_covariance_matrix(data, cov_matrix, m, dim);


    double *eigenvalues = (double *)malloc(dim * sizeof(double));
    double *eigenvectors = (double *)malloc(dim * dim * sizeof(double));
    eigen_decompose(cov_matrix, eigenvalues, eigenvectors, dim);


    double *projected_data = (double *)malloc(output_samples * k * sizeof(double));
    project_data(data, eigenvectors, projected_data, m, dim, k, sample_indices, output_samples);


    if (sample_indices == NULL) {

        for (int i = 0; i < output_samples; i++) {
            for (int j = 0; j < k; j++) {
                output[i][j] = projected_data[i * k + j];
            }
            output[i][k] = input[i][n - 1];
        }
    } else {

        for (int idx = 0; idx < output_samples; idx++) {
            int i = sample_indices[idx];
            for (int j = 0; j < k; j++) {
                output[idx][j] = projected_data[idx * k + j];
            }
            output[idx][k] = input[i][n - 1];
        }
    }


    free(data);
    free(cov_matrix);
    free(eigenvalues);
    free(eigenvectors);
    free(projected_data);
    if (sample_indices != NULL) {
        free(sample_indices);
    }
}


void PCA_transform(double **input, double **output, int m, int n, int k) {
    PCA_transform_advanced(input, output, m, n, k, -1, 0);
}


double** create_output_matrix(int rows, int cols) {
    double **matrix = (double**)malloc(rows * sizeof(double*));
    if (matrix == NULL) {
        return NULL;
    }
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*)calloc(cols, sizeof(double));
        if (matrix[i] == NULL) {

            for (int j = 0; j < i; j++) {
                free(matrix[j]);
            }
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}


void free_matrix(double **matrix, int rows) {
    if (matrix != NULL) {
        for (int i = 0; i < rows; i++) {
            if (matrix[i] != NULL) {
                free(matrix[i]);
            }
        }
        free(matrix);
    }
}