# Answer the following questions
1. If using the same pakage to express the special number, like pi, e, the integral is the true value.
2. It's the minimum error that cannot be reduced
3. The minimum error is 0
4. If using the numerical integration, it's close to the one in point 1. Relative error is 8.400431313389404e-07. While if using the package to integrate it, the result is the same as point 1.

# Sourve Code

## compiled language
C++ (GSL)
```C++
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_vector.h>
#include <vector>
#include <iostream>
#include <gsl/gsl_integration.h>
#include <iomanip>

using namespace std;

int N{10}; // number of points
double x_inf{0}; // inferior integration limit
double x_sup{M_PI/2}; // superior integration limit
double expected_result{(gsl_sf_exp(M_PI/2)-1)/2}; // expected result

int main (int argc, char *argv[]) {
    if (argc > 1 && argv[1]) N = atoi(argv[1]);
    if (argc > 2 && argv[2]) x_inf = atof(argv[2]);
    if (argc > 3 && argv[3]) x_sup = atof(argv[3]);

    double x;
    double range_x = x_sup - x_inf;
    double step_x = range_x / N;

    // Define the GSL function
    gsl_function* f = new gsl_function;
    f->function = [](double x, void* params) {
        return gsl_sf_cos(x)*gsl_sf_exp(x);
    };
    f->params = nullptr;

    // relative error tolerance
    auto epsrel_f = [](double I, double I_ex) {
        return fabs((I - I_ex) / I_ex);
    };

    gsl_vector* x_vector = gsl_vector_alloc(N+1);
    int i=0;
    // cout << "N= " << N << " x_inf= " << x_inf << " x_sup= " << x_sup << " step_x= " << step_x << endl;
    FILE* dp = fopen("./data_points.dat", "w");
    fprintf(dp, "x\t\t f(x)\n");
    for (x=x_inf; x<=x_sup; x+=step_x) {
        // cout << "x= " << x << " f(x)= " << f->function(x, f->params) << endl;
        double fx = f->function(x, f->params);
        gsl_vector_set(x_vector, i, fx);
        ++i;
        fprintf(dp, "%.4f \t %.4f\n", x, fx);
    }
    gsl_vector_free(x_vector);
    fclose(dp);

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
    double result, error;
    gsl_integration_qags(f, x_inf, x_sup, 0, 1e-7, 1000, w, &result, &error);

    cout << "Integration result: " << std::setprecision(16) << result << " ± " << error << endl;

    gsl_integration_workspace_free(w);
    delete f;

    cout << "Expected result: " << std::setprecision(16) << expected_result << endl;

    cout << "Relative error: " << std::setprecision(16) << epsrel_f(result, expected_result) << endl;

    return 0;
}
```
Makefile
```makefile
CXX = g++
CXXGSLFLAGS += `gsl-config --cflags` -std=c++17 -O2
LDGSLFLAGS  += `gsl-config --libs`

OBJ = function.o

TARGETS = function

function: function.cxx
	${CXX} ${CXXGSLFLAGS} -o $@ $< ${LDGSLFLAGS}

.PHONY: clean

clean:
	rm -f $(TARGETS)
```

## interpreted language
Python (SciPy)
```python
import numpy as np
from scipy.integrate import quad

expected_result = (np.exp(np.pi/2) - 1) / 2

def f(x):
    return np.exp(x) * np.cos(x)

def relative_error(approx, exact):
    return abs(approx - exact) / abs(exact)

if __name__ == "__main__":
    result, error = quad(f, 0, np.pi/2)
    print(f"Integral result: {result:.15f}, Estimated error: {error}")
    print(f"Expected result: {expected_result:.15f}")
    print(f"Relative error: {relative_error(result, expected_result)}")

    # integrate from .dat file
    with open("data_points.dat", "r") as file:
        data = np.loadtxt(file, skiprows=1)
    x_data, y_data = data[:, 0], data[:, 1]

    # Perform numerical integration using the trapezoidal rule
    result = np.trapezoid(y_data, x_data)
    print(f"Integral result from data: {result:.15f}")
    print(f"Relative error from data: {relative_error(result, expected_result)}")
```
