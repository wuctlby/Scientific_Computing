# Answer the questions
- [x] Verify the correctness of both results and benchmark the different solutions on both GPU and CPU platforms.
    - I simply followed userguide of CUDA to implement the DAXPY operation on GPU. The code is in [compiled language](#compiled-language) section below.
    ```
    Standard CPU DAXPY took 0.606833 seconds.
    Sum of result_cpu: -25044.4373296456
    OpenMP CPU DAXPY took 0.4113354270230047 seconds.
    Sum of result_cpu_omp: -25044.4373296456
    GPU DAXPY took 0.00544652795791626 seconds.
    Sum of result_gpu: -25044.4373296456
    Results match within tolerance.
    Results match within tolerance.
    ```
# Source Code

## compiled language
nvcc (cuda)
```nvcc
#include <iostream>
#include <vector>
#include <cmath>
#include <cuda_runtime.h>
#include <ctime>
#include <random>
#include <omp.h>
#include <iomanip>

using namespace std;

vector<double> daxpy_cpu(int n, double a, const vector<double>& x, const vector<double>& y) {
    vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = a * x[i] + y[i];
    }
    return result;
}

vector<double> daxpy_cpu_omp(int n, double a, const vector<double>& x, const vector<double>& y) {
    vector<double> result(n);
    #pragma omp parallel for schedule(static, 1024)
    for (int i = 0; i < n; i++) {
        result[i] = a * x[i] + y[i];
    }
    return result;
}

// https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#device-memory
__global__ void daxpy_kernel(int n, double a, const double* x, double* y, double* result_gpu) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        result_gpu[i] = a * x[i] + y[i];
    }
}

void verify_results(const double sum_a, const double sum_b) {
    double eps = numeric_limits<double>::epsilon() * 10;
    if (fabs(sum_a - sum_b) > eps) {
        cout << "Results do not match! Difference: " << fabs(sum_a - sum_b) << endl;
    } else {
        cout << "Results match within tolerance." << endl;
    }
}

int main() {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> normal_list(0.0L, 1.0L); // mean 0, stddev 1

    const int N = 100000000;
    const double a = 2.0L;
    vector<double> x(N), y(N);
    for (int i = 0; i < N; ++i) {
        x[i] = normal_list(gen);
        y[i] = normal_list(gen);
    }

    vector<double> result_cpu(N), result_cpu_omp(N), result_gpu(N);

    double start_time_cpu_standard = omp_get_wtime();
    result_cpu = daxpy_cpu(N, a, x, y);
    double end_time_cpu_standard = omp_get_wtime();
    cout << "Standard CPU DAXPY took " << (end_time_cpu_standard - start_time_cpu_standard) << " seconds." << endl;
    double sum_cpu = std::accumulate(result_cpu.begin(), result_cpu.end(), 0.0L);
    cout << "Sum of result_cpu: " << std::setprecision(16) << sum_cpu << endl;

    double start_time_cpu_omp = omp_get_wtime();
    result_cpu_omp = daxpy_cpu_omp(N, a, x, y);
    double end_time_cpu_omp = omp_get_wtime();
    cout << "OpenMP CPU DAXPY took " << (end_time_cpu_omp - start_time_cpu_omp) << " seconds." << endl;
    double sum_cpu_omp = std::accumulate(result_cpu_omp.begin(), result_cpu_omp.end(), 0.0L);
    cout << "Sum of result_cpu_omp: " << std::setprecision(16) << sum_cpu_omp << endl;

    // Allocate device memory and copy data from host to device
    size_t vector_size = N * sizeof(double);
    double *d_x, *d_y;
    cudaMalloc((void**)&d_x, vector_size);
    cudaMalloc((void**)&d_y, vector_size);
    cudaMemcpy(d_x, x.data(), vector_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, y.data(), vector_size, cudaMemcpyHostToDevice);
    double *d_result_gpu;
    cudaMalloc((void**)&d_result_gpu, vector_size);

    // Invoke GPU kernel
    int threads_per_block = 256;
    int blocks_per_grid = (N + threads_per_block - 1) / threads_per_block; // 为什么要加上 threads_per_block - 1？
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    daxpy_kernel<<<blocks_per_grid, threads_per_block>>>(N, a, d_x, d_y, d_result_gpu);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);

    cudaMemcpy(result_gpu.data(), d_result_gpu, vector_size, cudaMemcpyDeviceToHost);

    float runtime;
    cudaEventElapsedTime(&runtime, start, stop);
    cout << "GPU DAXPY took " << runtime / 1000.0 << " seconds." << endl;
    double sum_gpu = std::accumulate(result_gpu.begin(), result_gpu.end(), 0.0L);
    cout << "Sum of result_gpu: " << std::setprecision(16) << sum_gpu << endl;

    // Copy result back to host
    cudaMemcpy(result_gpu.data(), d_y, vector_size, cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(d_x);
    cudaFree(d_y);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    // Verify results
    verify_results(sum_cpu, sum_gpu);
    verify_results(sum_cpu_omp, sum_gpu);

    return 0;
}
```

Makefile
```makefile
NVCC = nvcc

CXXFLAGS  = -std=c++17 -O2 $(shell gsl-config --cflags)
LDFLAGS   = $(shell gsl-config --libs)
OMPFLAGS  = -Xcompiler -fopenmp 

TARGET = daxpy
SRC    = daxpy.cu

all: $(TARGET)

$(TARGET): $(SRC)
	$(NVCC) $(CXXFLAGS) $(OMPFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TARGET)
```