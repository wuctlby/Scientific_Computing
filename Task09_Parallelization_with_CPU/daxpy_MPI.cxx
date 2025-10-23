#include <iostream>
#include <vector>
#include <random>
#include <iomanip>
#include <cassert>
#include <time.h>
#include <omp.h>
#include <mpi.h>

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
    #pragma omp parallel for schedule(static, 1024)
    for (int i=0; i<N; ++i) {
        result[i] = a * x[i] + y[i];
    }
    return result;
}

vector<long double> my_daxpy_chunks_lb(const vector<long double>& y, const vector<long double>& x, double a, int chunk_size) {
    int N = y.size();
    vector<long double> result(N);

    int number_of_chunks = ceil(N / chunk_size);
    {
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

int main(int argc, char** argv) {
    // Initialize MPI first
    MPI_Init(&argc, &argv);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    auto result = vector<long double>();

    if (world_rank == 0) {
        cout << "Using standard random library from c++" << endl;
    }
    
    mt19937 r_generator(42);
    normal_distribution<double> normal_dist(0.0, 1.0); // mean 0, stddev 1
    vector<long double> x(N), y(N);
    for (int i=0; i<N; ++i) {
        x[i] = normal_dist(r_generator);
        y[i] = normal_dist(r_generator);
    }

    if (world_rank == 0) {
        double start_standard_time = omp_get_wtime();
        result = my_daxpy_lb(y, x, a);
        double end_standard_time = omp_get_wtime();
        double standard_duration = end_standard_time - start_standard_time;
        cout << "Standard DAXPY with long double took " << standard_duration << " seconds." <<  endl;
        cout << "Sum of standard result: " << std::setprecision(16) << std::accumulate(result.begin(), result.end(), 0.0L) << endl;

        double start_omp_time = omp_get_wtime();
        auto result_omp = my_daxpy_lb_omp(y, x, a);
        double end_omp_time = omp_get_wtime();
        double omp_duration = end_omp_time - start_omp_time;
        cout << "OpenMP DAXPY with long double took " << omp_duration << " seconds." <<  endl;
        cout << "Sum of result_omp: " << std::setprecision(16) << std::accumulate(result_omp.begin(), result_omp.end(), 0.0L) << endl;

        long double error_omp = 1;
        long double total_sum_omp = std::accumulate(result_omp.begin(), result_omp.end(), 0.0L);
        long double total_sum_standard = std::accumulate(result.begin(), result.end(), 0.0L);
        while (abs(total_sum_omp - total_sum_standard) < error_omp)
        {
            error_omp *= 0.1;
        }
        cout << "Final difference for OpenMP " << abs(total_sum_omp - total_sum_standard) << " is not less than epsilon " << error_omp << endl;
    }

    // MPI parallel computation
    double start_mpi_time = MPI_Wtime();
    
    // Calculate local chunk size for each process
    int local_N = N / world_size;
    int remainder = N % world_size;
    int local_start = world_rank * local_N + min(world_rank, remainder);
    if (world_rank < remainder) {
        local_N++;
    }
    
    // Each process works on its local portion
    vector<long double> local_result = my_daxpy_lb(
        vector<long double>(y.begin() + local_start, y.begin() + local_start + local_N),
        vector<long double>(x.begin() + local_start, x.begin() + local_start + local_N),
        a
    );
    
    double end_mpi_time = MPI_Wtime();
    double mpi_duration = end_mpi_time - start_mpi_time;
    
    // Gather results back to rank 0 (optional)
    vector<long double> result_mpi;
    if (world_rank == 0) {
        result_mpi.resize(N);
    }
    
    // Gather all local results to rank 0
    vector<int> recvcounts(world_size);
    vector<int> displs(world_size);
    MPI_Gather(&local_N, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (world_rank == 0) {
        displs[0] = 0;
        for (int i = 1; i < world_size; ++i) {
            displs[i] = displs[i-1] + recvcounts[i-1];
        }
    }
    
    MPI_Gatherv(local_result.data(), local_N, MPI_LONG_DOUBLE,
                result_mpi.data(), recvcounts.data(), displs.data(), 
                MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (world_rank == 0) {
        cout << "MPI DAXPY with long double (using " << world_size << " processes) took " 
             << mpi_duration << " seconds." << endl;
        cout << "Length of MPI result: " << result_mpi.size() << endl;
        cout << "Sum of MPI result: " << std::setprecision(16) << std::accumulate(result_mpi.begin(), result_mpi.end(), 0.0L) << endl;
    }

    if (world_rank == 0) {
        long double epsilon_mpi = 1;
        long double total_sum_standard = std::accumulate(result.begin(), result.end(), 0.0L);
        long double total_sum_mpi = std::accumulate(result_mpi.begin(), result_mpi.end(), 0.0L);
        while (abs(total_sum_mpi - total_sum_standard) < epsilon_mpi)
        {
            epsilon_mpi *= 0.1;
        }
        cout << "Final difference for MPI " << abs(total_sum_mpi - total_sum_standard) << " is not less than epsilon " << epsilon_mpi << endl;
    }
    
    MPI_Finalize();

    return 0;
}
