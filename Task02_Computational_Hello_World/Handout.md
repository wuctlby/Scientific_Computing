## Question:
### Did you find any problems in running the codes for some N. If so, do you have an idea why?
- the loop with large N in Python, if using the Standard list is quit slow. while using the `Numpy` is even faster than C++
- It's easy to make an mistake on the loop pf N*N matrix in C++, or forget to reset the matrix to get wrong results

### Where you able to test correctly the sum and product of points 1-3? If so, how? If not, what was the problem?
- Yes, I did. Using `np.all` in python to compare with the expected value. Comparing with the expected value one by one when calculating it. There might be an extended library in C++ that optimises matrix calculations, but I didn't use it.

## 02: Computational hello word

Implement computational "hello word" meme with interepted and compiled language

### Intereted language - Python

- Create an python enviroment and install dependencies by `uv`

    Step 1: create an enviroment
    ```bash
    python3 -m venv .venv # `.venv` must be used to satisfy the default env name of `uv`
    source .venv/bin/activate
    ```

    Step 2: install uv and needed dependencies
    ```bash
    pip install uv # install uv, the modden project manager tool
    uv init # initialize the project --> produce a `pyproject.toml`
    uv add numpy # or add 'numpy'  into [project.dependencies] part of `pyproject.toml`
    python -m pip show numpy # confirm we install the `numpy` successfully
    ```

- [source code](#intereted-language-python-python-source-code)

- Run `Task02_Computational_Hello_World/computational_hello_world.py` to perform the computational "hello word"
    ```bash
    # compare the performence of standard `Python List` and `Numpy` arrays
    python3 ./Task02_Computational_Hello_World/computational_hello_world.py -p 1
    ```

### Compiled language - C++

- Compile `Task02_Computational_Hello_World/computational_hello_world.cxx` to an executable `.o` file
    ```bash
    g++ -std=c++17 ./Task02_Computational_Hello_World/computational_hello_world.cxx -DCOMPILER_FLAGS="\"-std=c++17\"" -o computational_hello_world
    # `-std=c++17` use the C++ 17 standard libeliry
    ./computational_hello_world 1
    ```

- [source code](#compiled-language-c-c-source-code)

- Compile `Task02_Computational_Hello_World/computational_hello_world.cxx` with `-O2` (or `-O3`) to optimize the loop
    ```bash
    g++ -O2 -std=c++17 ./Task02_Computational_Hello_World/computational_hello_world.cxx -DCOMPILER_FLAGS="\"-O2 -std=c++17\"" -o computational_hello_world
    # `-O2` optimize the loop; `-std=c++17` use the C++ 17 standard libeliry
    ./computational_hello_world 1
    ```

### Matrix multiplication
- Interested language (only the one using numpy)
    ```bash
    python3 ./Task02_Computational_Hello_World/computational_hello_world.py -p 3
    ```
- Compiled language
    ```bash
    g++ -O2 -std=c++17 ./Task02_Computational_Hello_World/computational_hello_world.cxx -DCOMPILER_FLAGS="\"-O2 -std=c++17\"" -o computational_hello_world
    # `-O2` optimize the loop; `-std=c++17` use the C++ 17 standard libeliry
    ./computational_hello_world 3
    ```

### source code
#### Intereted language: Python {#python-source-code}
    ```python
    import argparse
    import time
    import numpy as np

    # scalar
    a = None

    # vector
    x = []
    y = []

    #results
    result_d = []

    # dimension
    N = [10, 10**6, 10**8]


    def compute_hello_world(a_value, x_value, y_value, dimension, equation, result_d_value=None):
        global a, x, y, result_d
        a = a_value
        x = [x_value] * dimension
        y = [y_value] * dimension

        for i in range(dimension):
            result_d.append(eval(equation, {'a': a, 'x': x[i], 'y': y[i]}))
            if result_d_value is not None and result_d[i] != result_d_value:
                raise ValueError(f"Computation error at index {i}: expected {result_d_value}, got {result_d[i]}")
        print(f"Computation completed for dimension {dimension}. All results match expected value: {result_d_value}")
        return result_d, a, x, y

    def p1():
        a_value = 3
        x_value = 0.1
        y_value = 7.1
        start_time = time.time()
        print("Starting computations with normal Python lists...")
        equation = 'a * x + y'
        result_d_value = eval(equation, {'a': a_value, 'x': x_value, 'y': y_value})
        with open("results_P1.txt", "w") as f:
            f.write(f"Equation: {equation}\n")
            f.write(f"Expected result for all elements: {result_d_value}\n")
        print(f"Expected result for equation '{equation}' for all elements: {result_d_value}")

        for dimension in N:
            results, a, x, y = compute_hello_world(a_value, x_value, y_value, dimension, equation, result_d_value)
            
            middle_time = time.time()
            print(f"Time taken for dimension {dimension}: {middle_time - start_time:.2f} seconds")
            with open("results_P1.txt", "a") as f:
                f.write(f"Dimension: {dimension}, All results match expected value: {result_d_value}\n")
                f.write(f"Time taken: {middle_time - start_time:.2f} seconds for dimension {dimension}\n")
            result_d = []

        end_time = time.time()
        print(f"All computations completed in {end_time - start_time:.2f} seconds \n")

        start_time = time.time()
        print("Starting computations with NumPy arrays...")
        a_value = 3
        x_value = 0.1
        y_value = 7.1
        equation = 'a * x + y'
        result_d_value = eval(equation, {'a': a_value, 'x': x_value, 'y': y_value})
        with open("results_P1.txt", "a") as f:
            f.write(f"\nUsing NumPy arrays:\n")
            f.write(f"Equation: {equation}\n")
            f.write(f"Expected result for all elements: {result_d_value}\n")
        print(f"Expected result for equation '{equation}' for all elements: {result_d_value}")
        for dimension in N:
            a = a_value
            x = np.full(dimension, x_value)
            y = np.full(dimension, y_value)
            result_d = eval(equation, {'a': a, 'x': x, 'y': y})
            if not np.all(result_d == result_d_value):
                raise ValueError(f"Computation error: not all results match expected value {result_d_value}")
            print(f"Computation completed for dimension {dimension}. All results match expected value: {result_d_value}")
            middle_time = time.time()
            print(f"Time taken for dimension {dimension}: {middle_time - start_time:.2f} seconds")
            with open("results_P1.txt", "a") as f:
                f.write(f"Dimension: {dimension}, All results match expected value: {result_d_value}\n")
                f.write(f"Time taken: {middle_time - start_time:.2f} seconds for dimension {dimension}\n")
            result_d = []
        end_time = time.time()
        print(f"All computations completed in {end_time - start_time:.2f} seconds")
        with open("results_P1.txt", "a") as f:
            f.write(f"All computations completed in {end_time - start_time:.2f} seconds\n")

    def p3():
        C = np.array([])
        A = np.array([])
        B = np.array([])
        N = [10, 100, 10000]

        A_value = 3
        B_value = 7.1

        start_time = time.time()
        expected_value = A_value * B_value
        print("Starting matrix multiplication computations with NumPy arrays...")
        print(f"Expected result for all elements: {expected_value}")
        with open("results_P3.txt", "w") as f:
            f.write(f"Matrix Multiplication A * B\n")
            f.write(f"Expected result for all elements: {expected_value}\n")
        for dimension in N:
            A = np.full((dimension, dimension), A_value)
            B = np.full((dimension, dimension), B_value)
            C = np.multiply(A, B)
            if not np.all(C == expected_value):
                raise ValueError(f"Computation error: not all results match expected value {expected_value}")
            print(f"Computation completed for dimension {dimension}. All results match expected value: {expected_value}")
            middle_time = time.time()
            print(f"Time taken for dimension {dimension}: {middle_time - start_time:.2f} seconds")
            with open("results_P3.txt", "a") as f:
                f.write(f"Dimension: {dimension}, All results match expected value: {expected_value}\n")
                f.write(f"Time taken: {middle_time - start_time:.2f} seconds for dimension {dimension}\n")
            C = []
        end_time = time.time()
        print(f"All computations completed in {end_time - start_time:.2f} seconds")
        with open("results_P3.txt", "a") as f:
            f.write(f"All computations completed in {end_time - start_time:.2f} seconds\n")

    if __name__ == "__main__":
        parser = argparse.ArgumentParser(description="Run computational hello world tasks.")
        parser.add_argument('--point', '-p', type=int, choices=[1, 3], required=True,
                            help="Specify which part to run: 1 for scalar/vector operations, 3 for matrix multiplication.")
        args = parser.parse_args()
        if args.point == 1:
            p1()
        elif args.point == 3:
            p3()
    ```

#### Compiled language: C++ {#C++-source-code}
    ```C++
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
    ```
