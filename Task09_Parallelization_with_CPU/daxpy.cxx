#include <iostream>
#include <vector>
#include <random>
#include <iomanip>
#include <cassert>
#include <time.h>
#include <omp.h>

using namespace std;

int N = 10000000; // size of vectors
int a = 2;  // scalar multiplier

vector<double> my_daxpy(const vector<double>& y, const vector<double>& x, double a) {
    int N = y.size();
    vector<double> result(N);
    for (int i=0; i<N; ++i) {
        result[i] = a * x[i] + y[i];
    }
    return result;
}

vector<long double> my_daxpy_lb(const vector<long double>& y, const vector<long double>& x, double a) {
    int N = y.size();
    vector<long double> result(N);
    for (int i=0; i<N; ++i) {
        result[i] = a * x[i] + y[i];
    }
    return result;
}

vector<long double> my_daxpy_lb_omp(const vector<long double>& y, const vector<long double>& x, double a) {
    int N = y.size();
    vector<long double> result(N);
    #pragma omp parallel for schedule(static)
    for (int i=0; i<N; ++i) {
        result[i] = a * x[i] + y[i];
    }
    return result;
}

vector<long double> my_daxpy_chunks_lb(const vector<long double>& y, const vector<long double>& x, double a, int chunk_size) {
    int N = y.size();
    vector<long double> result(N);

    int number_of_chunks = ceil(N / chunk_size);
    #pragma omp parallel
    {
        #pragma omp for
        for (int chunk_index = 0; chunk_index < number_of_chunks; ++chunk_index) {
            int current_index = chunk_index * chunk_size;
            int end_index = min(current_index + chunk_size, N);
            long double chunk_sum = 0.0;
            for (int i = current_index; i < end_index; ++i) {
                result[i] = a * x[i] + y[i];
                chunk_sum += result[i];
            }
        }
    }
    return result;
}

int main() {

    cout << "Using standard random library from c++" << endl;
    random_device rd;  // obtain a random number from hardware
    mt19937 r_generator(rd()); // seed the generator
    normal_distribution<double> normal_dist(0.0, 1.0); // mean 0, stddev 1
    vector<long double> x(N), y(N);
    for (int i=0; i<N; ++i) {
        x[i] = normal_dist(r_generator);
        y[i] = normal_dist(r_generator);
    }

    double start_standard_time = omp_get_wtime();
    auto result = my_daxpy_lb(y, x, a);
    double end_standard_time = omp_get_wtime();
    double standard_duration = end_standard_time - start_standard_time;
    cout << "Standard DAXPY with long double took " << standard_duration << " seconds." <<  endl;

    double start_chunked_time = omp_get_wtime();
    auto result_chunks = my_daxpy_lb_omp(y, x, a);
    double end_chunked_time = omp_get_wtime();
    double chunked_duration = end_chunked_time - start_chunked_time;
    cout << "OpenMP DAXPY with long double took " << chunked_duration << " seconds." <<  endl;

    return 0;
}
