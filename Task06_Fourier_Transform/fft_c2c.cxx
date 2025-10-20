#include <random>
#include <complex>
#include <vector>
#include <iostream>
#include <fftw3.h>
#include <algorithm>

using namespace std;

/// @brief Save matrix to file
/// @tparam T data type (fftw_complex or double)
/// @param filename name of the file
/// @param matrix pointer to the matrix data
/// @param N size of the matrix (N x N)
/// @param scale scaling factor for the output values
template<typename T>
void save_matrix_to_file(const char* filename, T* matrix, int N, double scale = 1.0, bool doR2C = false) {
    FILE* fp = fopen(filename, "w");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int idx = i * N + j;
            if constexpr (is_same<T, fftw_complex>::value) {
                if (doR2C) {
                    if (j >= N / 2 + 1) {
                        fprintf(fp, "% 12.6f + % 12.6fi\t\t", 0.0, 0.0); // zero padding for real-to-complex FFT output
                    } else {
                        fprintf(fp, "% 12.6f + % 12.6fi\t\t", matrix[idx - (i * (N / 2 - 1))][0] * scale, matrix[idx - (i * (N / 2 - 1))][1] * scale);
                    }
                } else {
                    fprintf(fp, "% 12.6f + % 12.6fi\t\t", matrix[idx][0] * scale, matrix[idx][1] * scale);
                }
            } else if constexpr (is_same<T, double>::value) {
                fprintf(fp, "% 12.6f\t\t", matrix[idx] * scale);
            }
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

// void save_matrix_to_file(const char* filename, fftw_complex* matrix, int N, double scale = 1.0) {
//     FILE* fp = fopen(filename, "w");
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < N; j++) {
//             int idx = i * N + j;
//             fprintf(fp, "% 12.6f + % 12.6fi\t\t", matrix[idx][0] * scale, matrix[idx][1] * scale);
//         }
//         fprintf(fp, "\n");
//     }
//     fclose(fp);
// }

template<typename T>
void mean_med_abs_rel_err(T* original, T* reconstructed, int N) {
    double mean_abs_err = 0.0;
    double mean_rel_err = 0.0;
    double med_abs_err = 0.0;
    double med_rel_err = 0.0;
    vector<double> abs_errors;
    vector<double> rel_errors;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int idx = i * N + j;
            if constexpr (is_same<T, fftw_complex>::value) {
                double orig_real = original[idx][0];
                double orig_imag = original[idx][1];
                double orig_magnitude = sqrt(pow(orig_real, 2) + pow(orig_imag, 2));
                double rec_real = reconstructed[idx][0] / (N * N);
                double rec_imag = reconstructed[idx][1] / (N * N);
                double abs_err = sqrt(pow(orig_real - rec_real, 2) + pow(orig_imag - rec_imag, 2));
                double rel_err = (orig_magnitude != 0.0) ? (abs_err / orig_magnitude) : 0.0;
                mean_abs_err += abs_err;
                mean_rel_err += rel_err;
                abs_errors.push_back(abs_err);
                rel_errors.push_back(rel_err);
            } else {
                double orig_value = original[idx];
                double rec_value = reconstructed[idx] / (N * N);
                double abs_err = fabs(orig_value - rec_value);
                double rel_err = (orig_value != 0.0) ? (abs_err / fabs(orig_value)) : 0.0;
                mean_abs_err += abs_err;
                mean_rel_err += rel_err;
                abs_errors.push_back(abs_err);
                rel_errors.push_back(rel_err);
            }
        }
    }
    mean_abs_err /= (N * N);
    mean_rel_err /= (N * N);
    nth_element(abs_errors.begin(), abs_errors.begin() + (N * N) / 2, abs_errors.end());
    med_abs_err = abs_errors[(N * N) / 2];
    nth_element(rel_errors.begin(), rel_errors.begin() + (N * N) / 2, rel_errors.end());
    med_rel_err = rel_errors[(N * N) / 2];
    cout << "Mean Absolute Error: " << mean_abs_err << endl;
    cout << "Median Absolute Error: " << med_abs_err << endl;
    cout << "Mean Relative Error: " << mean_rel_err << endl;
    cout << "Median Relative Error: " << med_rel_err << endl;
}

// void mean_med_abs_rel_err(fftw_complex* original, fftw_complex* reconstructed, int N) {
//     double mean_abs_err = 0.0;
//     double mean_rel_err = 0.0;
//     double med_abs_err = 0.0;
//     double med_rel_err = 0.0;
//     vector<vector<double>> errors;
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < N; j++) {
//             int idx = i * N + j;
//             double orig_real = original[idx][0];
//             double orig_imag = original[idx][1];
//             double orig_magnitude = sqrt(pow(orig_real, 2) + pow(orig_imag, 2));
//             double rec_real = reconstructed[idx][0] / (N * N);
//             double rec_imag = reconstructed[idx][1] / (N * N);

//             double abs_err = sqrt(pow(orig_real - rec_real, 2) + pow(orig_imag - rec_imag, 2));
//             double rel_err = (orig_magnitude != 0.0) ? (abs_err / orig_magnitude) : 0.0;
//             mean_abs_err += abs_err;
//             mean_rel_err += rel_err;
//             errors.push_back({abs_err, rel_err});
//         }
//     }
//     mean_abs_err /= (N * N);
//     mean_rel_err /= (N * N);
//     nth_element(errors[0].begin(), errors[0].begin() + N / 2, errors[0].end());
//     med_abs_err = errors[0][N / 2];
//     nth_element(errors[1].begin(), errors[1].begin() + N / 2, errors[1].end());
//     med_rel_err = errors[1][N / 2];
//     cout << "Mean Absolute Error: " << mean_abs_err << endl;
//     cout << "Median Absolute Error: " << med_abs_err << endl;
//     cout << "Mean Relative Error: " << mean_rel_err << endl;
//     cout << "Median Relative Error: " << med_rel_err << endl;

// }

int main() {
    // random number generator
    random_device rd;
    mt19937 r_generator(rd()); // seed the generator
    normal_distribution<double> normal_dist(0.0, 1.0); // mean 0, stddev 1

    int N = 6; // size of the array
    fftw_complex *matrix_A, *matrix_C, *matrix_inverseA;
    matrix_A = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    matrix_inverseA = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    matrix_C = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);

    cout << "Filling matrix with random complex numbers..." << endl;
    // Fill matrix_A with random complex numbers
    double sum_real = 0.0;
    double sum_imag = 0.0;
    for (int i = 0; i < N * N; ++i) {
        matrix_A[i][0] = normal_dist(r_generator); // real part
        matrix_A[i][1] = 0.; // imaginary part
        sum_real += matrix_A[i][0];
        sum_imag += matrix_A[i][1];
    }
    cout << "Sum of real parts: " << sum_real << endl;
    cout << "Sum of imaginary parts: " << sum_imag << endl;
    cout << "Average value of real part: " << sum_real / (N * N) << endl;
    cout << "Average value of imaginary part: " << sum_imag / (N * N) << endl;

    save_matrix_to_file<fftw_complex>("input_matrix.txt", matrix_A, N);

    // Create FFTW plan for forward transform
    fftw_plan forward;
    forward = fftw_plan_dft_2d(N, N, matrix_A, matrix_C, FFTW_FORWARD, FFTW_ESTIMATE);

    // Execute the FFT
    fftw_execute(forward);

    save_matrix_to_file<fftw_complex>("output_matrix.txt", matrix_C, N);

    fftw_plan backward;
    backward = fftw_plan_dft_2d(N, N, matrix_C, matrix_inverseA, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Execute the inverse FFT
    fftw_execute(backward);

    save_matrix_to_file<fftw_complex>("inverse_output_matrix.txt", matrix_inverseA, N, 1.0 / (N * N));

    // Calculate and print mean and median absolute and relative errors
    mean_med_abs_rel_err<fftw_complex>(matrix_A, matrix_inverseA, N);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    cout << "Performing real-to-complex and complex-to-real FFT..." << endl;
    // type conversion: complex to real
    double* matrix_realA = (double*) fftw_alloc_real(N * N);
    for (int i = 0; i < N * N; ++i) {
        matrix_realA[i] = matrix_A[i][0]; // take only the real part
    }
    save_matrix_to_file<double>("input_matrix_r2c.txt", matrix_realA, N);

    fftw_complex* matrix_R = (fftw_complex*) fftw_alloc_complex(N * (N / 2 + 1));
    double* matrix_inverseR = (double*) fftw_alloc_real(N * N);
    fftw_plan r2c_forward, c2r_backward;
    
    r2c_forward = fftw_plan_dft_r2c_2d(N, N, matrix_realA, matrix_R, FFTW_ESTIMATE);
    fftw_execute(r2c_forward);
    save_matrix_to_file<fftw_complex>("output_matrix_r2c.txt", matrix_R, N, 1.0, true);

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // obtain matrix C from matrix R
    vector<fftw_complex> matrix_C_from_R(N * N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N / 2 + 1; ++j) {
            int idx_R = i * (N / 2 + 1) + j;
            int idx_C = i * N + j;
            matrix_C_from_R[idx_C][0] = matrix_R[idx_R][0];
            matrix_C_from_R[idx_C][1] = matrix_R[idx_R][1];
            cout << "i: " << i << ", j: " << j << ", idx_C: " << idx_C << ", idx_R: " << idx_R << endl;
            if (j != 0 && j != N / 2) { // fill the conjugate symmetric part
                int idx_C_conj = i * N + (N - j);
                cout << "i: " << i << ", j: " << j << ", idx_C_conj: " << idx_C_conj << ", idx_R: " << idx_R << endl;
                matrix_C_from_R[idx_C_conj][0] = matrix_R[idx_R][0];
                matrix_C_from_R[idx_C_conj][1] = -matrix_R[idx_R][1];
            }
        }
    }
    save_matrix_to_file<fftw_complex>("output_matrix_r2c_conj.txt", matrix_C_from_R.data(), N);
    /////////////////////////////////////////////////////////////////////////////////////////////////

    c2r_backward = fftw_plan_dft_c2r_2d(N, N, matrix_R, matrix_inverseR, FFTW_ESTIMATE);
    fftw_execute(c2r_backward);
    save_matrix_to_file<double>("inverse_output_matrix_r2c.txt", matrix_inverseR, N, 1.0 / (N * N));

    mean_med_abs_rel_err<double>(matrix_realA, matrix_inverseR, N);

    // Cleanup
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
    fftw_destroy_plan(r2c_forward);
    fftw_destroy_plan(c2r_backward);
    fftw_free(matrix_A);
    fftw_free(matrix_C);
    fftw_free(matrix_inverseA);
    fftw_free(matrix_R);
    fftw_free(matrix_inverseR);
    fftw_free(matrix_realA);

    fftw_cleanup();

    cout << "Machine precision: " << endl;
    cout << "float epsilon: " << numeric_limits<float>::epsilon() << endl;
    cout << "double epsilon: " << numeric_limits<double>::epsilon() << endl;
    cout << "long double epsilon: " << numeric_limits<long double>::epsilon() << endl;
    return 0;
}