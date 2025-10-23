# Answer the questions
- [x] Try implementing some unit test for the code
    - 5 test cases added, see the code below
- [x] (Optional) do it for both interpreted and compiled language
  - see output below

# Source Code

## compiled language
C++ (GSL)
```C++
#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
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

void TEST_MY_DAXPY() {
    // Test case 1
    {
        vector<double> x = {1.0, 2.0, 3.0};
        vector<double> y = {4.0, 5.0, 6.0};
        double a = 2.0;
        vector<double> expected = {6.0, 9.0, 12.0};
        vector<double> result = my_daxpy(y, x, a);
        assert(result == expected);
    }
    // Test case 2
    {
        vector<double> x = {0.0, -1.0, 1.0};
        vector<double> y = {1.0, 1.0, 1.0};
        double a = 3.0;
        vector<double> expected = {1.0, -2.0, 4.0};
        vector<double> result = my_daxpy(y, x, a);
        assert(result == expected);
    }
    // Test case 3
    {
        vector<double> x = {1.5, 2.5};
        vector<double> y = {3.5, 4.5};
        double a = 0.0;
        vector<double> expected = {3.5, 4.5};
        vector<double> result = my_daxpy(y, x, a);
        assert(result == expected);
    }
    // Test case 4
    {
        vector<double> x = {};
        vector<double> y = {};
        double a = 5.0;
        vector<double> expected = {};
        vector<double> result = my_daxpy(y, x, a);
        assert(result == expected);
    }
    // Test case 5
    {
        vector<double> x = {1.0};
        vector<double> y = {2.0};
        double a = -1.0;
        vector<double> expected = {1.0};
        vector<double> result = my_daxpy(y, x, a);
        assert(result == expected);
    }
    cout << "All tests passed!" << endl;
}

int main() {
    TEST_MY_DAXPY();
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




    unsigned int r_generator_seed = rd(); // seed for GSL random generator
    cout << "Using GSL random library" << endl;
    const gsl_rng_type* T;
    gsl_rng* g_generator;
    T = gsl_rng_default;
    g_generator = gsl_rng_alloc(T);
    gsl_rng_set(g_generator, r_generator_seed);
    // seed the generator using current time and draw one Gaussian sample (sigma = 1.0)
    double gsl_sample = gsl_ran_gaussian(g_generator, 1.0);
    for (int i=0; i<N; ++i) {
        x[i] = gsl_ran_gaussian(g_generator, 1.0);
        y[i] = gsl_ran_gaussian(g_generator, 1.0);
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

    result = my_daxpy(y, x, a);
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

    gsl_rng_free(g_generator);

}
```

Makefile
```makefile
CXX = g++
CXXGSLFLAGS += `gsl-config --cflags` -std=c++17 -O2
LDGSLFLAGS  += `gsl-config --libs`

TARGETS = $(patsubst %.cxx, %, $(wildcard *.cxx))

All: $(TARGETS)

%: %.cxx
	${CXX} ${CXXGSLFLAGS} -o $@ $< ${LDGSLFLAGS}
```

## interpreted language
Python (NumPy), command: `python daxpy.py` and `pytest daxpy.py`
```python
import numpy as np
import pytest

def test_daxpy():
    # case 1
    a = 2.0
    x = np.array([1.0, 2.0, 3.0])
    y = np.array([4.0, 5.0, 6.0])
    expected = np.array([6.0, 9.0, 12.0])
    result = daxpy(a, x, y)
    np.testing.assert_allclose(result, expected)
    assert result.all() == expected.all()

    # case 2
    a = 3.0
    x = np.array([0.0, -1.0, 1.0])
    y = np.array([1.0, 1.0, 1.0])
    expected = np.array([1.0, -2.0, 4.0])
    result = daxpy(a, x, y)
    np.testing.assert_allclose(result, expected)
    assert result.all() == expected.all()

    # case 3
    a = 0.0
    x = np.array([1.5, 2.5])
    y = np.array([3.5, 4.5])
    expected = np.array([3.5, 4.5])
    result = daxpy(a, x, y)
    np.testing.assert_allclose(result, expected)
    assert result.all() == expected.all()

    # case 4
    a = 5.0
    x = np.array([])
    y = np.array([])
    expected = np.array([])
    result = daxpy(a, x, y)
    np.testing.assert_allclose(result, expected)
    assert result.all() == expected.all()
    print("All test cases passed!")

    # case 5 wrong onput
    a = -1.0
    x = np.array([1.0])
    y = np.array([2.0])
    expected = np.array([2.0])
    result = daxpy(a, x, y)
    np.testing.assert_allclose(result, expected)
    assert result.all() == expected.all()


def daxpy(a, x, y):
    if x.shape != y.shape:
        raise ValueError("Input vectors must have the same shape.")
    return a * x + y

if __name__ == "__main__":
    N = 1000
    # seed generator from device
    rng = np.random.default_rng()
    # generate random numbers from standard gaussian distribution
    x = rng.standard_normal((0, 1, N))
    y = rng.standard_normal((0, 1, N))
    a = rng.uniform(-10, 10)
    result = daxpy(a, x, y)
    print("DAXPY computation completed.")
    # Run tests
    test_daxpy()
```