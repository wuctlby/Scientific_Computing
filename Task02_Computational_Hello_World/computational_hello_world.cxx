#include <iostream>
#include <vector>
#include <time.h>
#include <functional>

using std::cout;
using std::endl;
using std::vector;
using std::function;
using std::string;
using std::stoi;

// scalar
double a;
// vectors
vector<double> x;
vector<double> y;
vector<double> result_d;
 // dimensions
 vector<int> N = {10, 1000000, 100000000};

void p2() {
    a = 3.0;
    double x_value = 0.1;
    double y_value = 7.1;

    FILE* results_txt = fopen("results_P2.txt", "a");
    if (results_txt == NULL) {
        cout << "Error opening file!" << endl;
    }
    clock_t start_time = clock();
    cout << "Hello World from C++!" << endl;
    fprintf(results_txt, "Hello World from C++!\n");

    #ifdef COMPILER_FLAGS
        cout << "Compiler flags: " << COMPILER_FLAGS << endl;
        fprintf(results_txt, "Compiler flags: %s\n", COMPILER_FLAGS);
    #else
        cout << "No compiler flags provided." << endl;
        fprintf(results_txt, "No compiler flags provided.\n");
    #endif

    function<double(double, double, double)> equation = [](double a, double x, double y) { return a * x + y; };

    double expected = equation(a, x_value, y_value);
    cout << "Expected result: " << expected << endl;
    fprintf(results_txt, "Expected result: %f\n", expected);

    for (int n : N) {
        x.resize(n, x_value);
        y.resize(n, y_value);
        result_d.resize(n);

        clock_t loop_start_time = clock();
        for (int i = 0; i < n; ++i) {
            result_d[i] = equation(a, x[i], y[i]);
            if (result_d[i] != expected) {
                cout << "Error at index " << i << ": got " << result_d[i] << ", expected " << expected << endl;
                fprintf(results_txt, "Error at index %d: got %f, expected %f\n", i, result_d[i], expected);
            }
        }
        clock_t loop_end_time = clock();
        double elapsed_time = double(loop_end_time - loop_start_time) / CLOCKS_PER_SEC;
        cout << "Time taken for n=" << n << ": " << elapsed_time << " seconds" << endl;
        fprintf(results_txt, "Time taken for n=%d: %f seconds\n", n, elapsed_time);
    }
    clock_t end_time = clock();
    double total_elapsed_time = double(end_time - start_time) / CLOCKS_PER_SEC;
    cout << "Total time taken: " << total_elapsed_time << " seconds" << endl;
    fprintf(results_txt, "Total time taken: %f seconds\n", total_elapsed_time);
    fclose(results_txt);

}

void p3() {
    vector<vector<double>> A;
    vector<vector<double>> B;
    vector<vector<double>> C;
    vector<int> N = {10, 100, 10000};
    double A_value = 3.0;
    double B_value = 7.1;
    double expected_value = A_value * B_value;
    cout << "Starting matrix multiplication computations..." << endl;
    cout << "Expected result for all elements: " << expected_value << endl;
    FILE* results_txt = fopen("results_P3.txt", "a");
    fwprintf(results_txt, L"Matrix Multiplication A * B from C++ ---\n");
    fwprintf(results_txt, L"Expected result for all elements: %f\n", expected_value);
    clock_t start_time = clock();
    for (int dimension : N) {
        A.resize(dimension, vector<double>(dimension, A_value));
        B.resize(dimension, vector<double>(dimension, B_value));
        C.resize(dimension, vector<double>(dimension, 0.0));
        clock_t loop_start_time = clock();
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                C[i][j] = A[i][j] * B[j][i];
                if (C[i][j] != expected_value) {
                    cout << "Error at position (" << i << ", " << j << "): got " << C[i][j] << ", expected " << expected_value << endl;
                    fwprintf(results_txt, L"Error at position (%d, %d): got %f, expected %f\n", i, j, C[i][j], expected_value);
                }
            }
        }
        clock_t loop_end_time = clock();
        double elapsed_time = double(loop_end_time - loop_start_time) / CLOCKS_PER_SEC;
        cout << "Time taken for dimension=" << dimension << ": " << elapsed_time << " seconds" << endl;
        fwprintf(results_txt, L"Time taken for dimension=%d: %f seconds\n", dimension, elapsed_time);
        A.clear();
        B.clear();
        C.clear();
    }
    clock_t end_time = clock();
    double total_elapsed_time = double(end_time - start_time) / CLOCKS_PER_SEC;
    cout << "Total time taken: " << total_elapsed_time << " seconds" << endl;
    fwprintf(results_txt, L"Total time taken: %f seconds\n", total_elapsed_time);
    fclose(results_txt);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: ./program <point>" << endl;
        return 1;
    }

    int point = stoi(argv[1]);

    if (point == 2) {
        p2();
    } else if (point == 3) {
        p3();
    } else {
        cout << "Invalid point specified. Use 2 or 3." << endl;
    }

    return 0;
}
