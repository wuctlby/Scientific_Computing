#include <iostream>
#include <vector>
#include <random>
#include <iomanip>
#include <cassert>

using namespace std;

int N = 100000; // size of vectors
int a = 2;  // scalar multiplier

vector<double> my_daxpy(const vector<double>& y, const vector<double>& x, double a) {
    int N = y.size();
    vector<double> result(N);
    for (int i=0; i<N; ++i) {
        result[i] = a * x[i] + y[i];
    }
    return result;
}

pair<vector<double>, vector<double>> my_daxpy_chunks(const vector<double>& y, const vector<double>& x, double a, int chunk_size) {
    int N = y.size();
    vector<double> result(N);
    vector<double> partial_chunk_sums;
    int number_of_chunks = ceil(N / chunk_size);
    for (int chunk_index = 0; chunk_index < number_of_chunks; ++chunk_index) {
        int current_index = chunk_index * chunk_size;
        int end_index = min(current_index + chunk_size, N);
        double chunk_sum = 0.0;
        for (int i = current_index; i < end_index; ++i) {
            result[i] = a * x[i] + y[i];
            chunk_sum += result[i];
        }
        partial_chunk_sums.push_back(chunk_sum);
    }
    return make_pair(result, partial_chunk_sums);
}

pair<vector<long double>, vector<long double>> my_daxpy_chunks_lb(const vector<double>& y, const vector<double>& x, double a, int chunk_size) {
    int N = y.size();
    vector<long double> result(N);
    vector<long double> partial_chunk_sums;

    int number_of_chunks = ceil(N / chunk_size);
    for (int chunk_index = 0; chunk_index < number_of_chunks; ++chunk_index) {
        int current_index = chunk_index * chunk_size;
        int end_index = min(current_index + chunk_size, N);
        long double chunk_sum = 0.0;
        for (int i = current_index; i < end_index; ++i) {
            result[i] = a * x[i] + y[i];
            chunk_sum += result[i];
        }
        partial_chunk_sums.push_back(chunk_sum);
    }
    return make_pair(result, partial_chunk_sums);
}

int main() {
    cout << "Using standard random library from c++" << endl;
    random_device rd;  // obtain a random number from hardware
    mt19937 r_generator(rd()); // seed the generator
    normal_distribution<double> normal_dist(0.0, 1.0); // mean 0, stddev 1
    vector<double> x(N), y(N);
    for (int i=0; i<N; ++i) {
        x[i] = normal_dist(r_generator);
        y[i] = normal_dist(r_generator);
    }

    if (N <= 15) {
        cout << "Vector x: \n";
        for (const auto& val : x) cout << fixed << setprecision(6) << val << " ";
        cout << endl;
        cout << "Mean of x: " << std::setprecision(16) << std::accumulate(x.begin(), x.end(), 0.0) / N << endl;
        cout << "Sum of x: " << std::setprecision(16) << std::accumulate(x.begin(), x.end(), 0.0) << endl;
        cout << "Vector y: \n";
        for (const auto& val : y) cout << fixed << setprecision(6) << val << " ";
        cout << endl;
        cout << "Mean of y: " << std::setprecision(16) << std::accumulate(y.begin(), y.end(), 0.0) / N << endl;
        cout << "Sum of y: " << std::setprecision(16) << std::accumulate(y.begin(), y.end(), 0.0) << endl;
    } else {
        cout << "Mean of x: " << std::setprecision(16) << std::accumulate(x.begin(), x.end(), 0.0) / N << endl;
        cout << "Sum of x: " << std::setprecision(16) << std::accumulate(x.begin(), x.end(), 0.0) << endl;
        cout << "Mean of y: " << std::setprecision(16) << std::accumulate(y.begin(), y.end(), 0.0) / N << endl;
        cout << "Sum of y: " << std::setprecision(16) << std::accumulate(y.begin(), y.end(), 0.0) << endl;
    }

    vector<double> result = my_daxpy(y, x, a);
    if (N <= 15) {
        cout << "Result of DAXPY (a*x + y) with a = " << a << ": \n";
        for (const auto& val : result) cout << fixed << setprecision(6) << val << " ";
        cout << endl;
        cout << "expected result (a*x + y): \n";
        for (int i=0; i<N; ++i) cout << fixed << setprecision(6) << a*x[i] + y[i] << " ";
        cout << endl;
    }
    cout << "Mean of result: " << std::setprecision(16) << std::accumulate(result.begin(), result.end(), 0.0) / N << endl;
    cout << "Sum of result: " << std::setprecision(16) << std::accumulate(result.begin(), result.end(), 0.0) << endl;
    cout << string(45, '-') << "\n" << endl;

    int chunk_size = 8;
    auto [result_chunks, partial_sums] = my_daxpy_chunks(y, x, a, chunk_size);
    assert(result == result_chunks);
    if (result_chunks == result) {
        cout << "Chunked DAXPY result matches the standard DAXPY result." << endl;
    } else {
        cout << "Chunked DAXPY result does NOT match the standard DAXPY result!" << endl;
    }

    long double error = 0.0;
    for (int i = 0; i < N; ++i) {
        error += abs(result[i] - result_chunks[i]);
    }
    cout << "Total absolute error between standard and chunked DAXPY results: " << error << endl;

    double total_sum_from_chunks = accumulate(partial_sums.begin(), partial_sums.end(), 0.0);
    cout << "Total sum from chunked DAXPY: " << total_sum_from_chunks << endl;
    
    double total_sum_from_chunks_result = accumulate(result_chunks.begin(), result_chunks.end(), 0.0);
    cout << "Total sum from chunked DAXPY result vector: " << total_sum_from_chunks_result << endl;

    double total_sum_standard = std::accumulate(result.begin(), result.end(), 0.0);
    cout << "Total sum from standard DAXPY: " << total_sum_standard << endl;
    // assert(total_sum_from_chunks == total_sum_standard);
    long double epsilon = 1;
    while (abs(total_sum_from_chunks - total_sum_standard) < epsilon)
    {
        epsilon *= 0.1;
    }
    cout << "Final difference " << abs(total_sum_from_chunks - total_sum_standard) << " is not less than epsilon " << epsilon << endl;
    long double epsilon_result = 1;
    while (abs(total_sum_from_chunks_result - total_sum_standard) < epsilon)
    {
        epsilon *= 0.1;
    }
    cout << "Final difference from result vector " << abs(total_sum_from_chunks_result - total_sum_standard) << " is not less than epsilon " << epsilon << endl;

    auto [result_chunks_lb, partial_sums_lb] = my_daxpy_chunks_lb(y, x, a, chunk_size);
    long double total_sum_from_chunks_lb = accumulate(partial_sums_lb.begin(), partial_sums_lb.end(), 0.0L);
    cout << "Total sum from chunked DAXPY with long double: " << total_sum_from_chunks_lb << endl;
    if (abs(total_sum_from_chunks_lb - total_sum_standard) < 1e-10) {
        cout << "Total sums with long double match!" << endl;
    } else {
        cout << "Total sums with long double do NOT match!" << endl;
    }
    long double difference_lb = abs(total_sum_from_chunks_lb - total_sum_standard);
    long double difference = abs(total_sum_from_chunks - total_sum_standard);
    long double difference_chunks = abs(total_sum_from_chunks_lb - total_sum_from_chunks);
    cout << "Difference between chunked long double and chunked double: " << difference_chunks << endl;
    cout << "Difference with long double: " << difference_lb << endl;
    cout << "Difference with double: " << difference << endl;

}
