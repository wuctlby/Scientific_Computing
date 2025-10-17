# Answer the questions
1. Not the same. I guess it's because of the floating point precision.
2. Testing it according to the mean value of the result vector whether close to 0.

# Sourve Code

## compiled language
C++ (GSL)
- Question 1
```C++
#include <vector>
#include <gsl/gsl_vector.h>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

vector<double> values = {1.0, 1.0e16, -1.0e16, -0.5};
vector<long double> values_ld = {1.0, 1.0e16, -1.0e16, -0.5};
int N = values.size();
int main () {
    cout << "Using loop to sum values: " << endl;
    cout << string(45, '-') << endl;
    double sum = 0.0;
    long double sum_ld = 0.0;
    for (int i=0; i<N; ++i) {
        sum += values[i];
        sum_ld += values_ld[i];
        cout << "After adding " << values[i] << ", " << "sum = " << sum << endl; // << fixed << setprecision(32) << "sum =  " << sum << endl;
        cout << "After adding " << values_ld[i] << ", " << "sum_ld = " << sum_ld << endl; // << fixed << setprecision(32) << "sum_ld =  " << sum_ld << endl;
    }
    cout << "Final sum = " << sum << endl; // << fixed << setprecision(32) << sum << endl;
    cout << "Final sum_ld = " << sum_ld << endl; // << fixed << setprecision(32) << sum_ld << endl;
    cout << string(45, '-') << "\n" << endl;

    cout << "Using GSL vector to sum values: " << endl;
    cout << string(45, '-') << endl;
    gsl_vector* v = gsl_vector_alloc(N);
    for (int i=0; i<N; ++i) {
        gsl_vector_set(v, i, values[i]);
    }
    double sum_gsl = gsl_vector_sum(v);
    cout << "sum_gsl = " << sum_gsl << endl; // << fixed << setprecision(32) << sum_gsl << endl;
    gsl_vector_free(v);
    cout << string(45, '-') << endl;
    cout << "\n" << endl;

    cout << "Using Kahan summation to sum values: " << endl;
    cout << string(45, '-') << endl;
    double sum_kahan = 0.0;
    double c = 0.0; // for lost low-order bits
    double t = 0.0; // sum
    double y = 0.0; // value with low-order bits corrected
    for (int i=0; i<N; ++i) {
        y = values[i] - c;    // consider the lost low-order bits
        cout << "Adding " << fixed << setprecision(32) << values[i] << ", y = " << fixed << setprecision(32) << y << " with low-order bits c = " << fixed << setprecision(32) << c << endl;
        t = sum_kahan + y;    // add y to sum
        cout << "Intermediate sum t = " << fixed << setprecision(32) << t << endl;
        if ( t == sum_kahan ) {
            cout << "Warning!!!!!!!!!!!!!! "<< endl;
            cout << "the precision of sum_kahan is missing the precision of y" << endl;
        } else {
            cout << "No precision is lost." << endl;
            cout << "But t = " << fixed << setprecision(32) << t << " sum_kahan = " << fixed << setprecision(32) << sum_kahan << " y = " << fixed << setprecision(32) << y << endl;
        }
        c = (t - sum_kahan) - y;    // record the lost low-order bits
        cout << "New low-order bits c = " << fixed << setprecision(32) << c << endl;
        sum_kahan = t;                // update sum
        cout << "Updated sum_kahan = " << fixed << setprecision(32) << sum_kahan << endl;
        cout << "---------------------------------------------" << endl;
    }
    cout << "Final sum_kahan = " << fixed << setprecision(32) << sum_kahan << " sum = " << fixed << setprecision(32) << t << endl;
    cout << string(45, '-') << endl;
    cout << "\n" << endl;

    cout << "Using long double Kahan summation to sum values: " << endl;
    cout << string(45, '-') << endl;
    long double sum_kahan_ld = 0.0;
    long double c_ld = 0.0; // for lost low-order bits
    long double t_ld = 0.0; // sum
    long double y_ld = 0.0; // value with low-order bits corrected
    for (int i=0; i<N; ++i) {
        y_ld = values_ld[i] - c_ld;    // consider the lost low-order bits
        cout << "Adding " << fixed << setprecision(32) << values_ld[i] << ", y_ld = " << fixed << setprecision(32) << y_ld << " with low-order bits c_ld = " << fixed << setprecision(32) << c_ld << endl;
        t_ld = sum_kahan_ld + y_ld;    // add y to sum
        cout << "Intermediate sum t_ld = " << fixed << setprecision(32) << t_ld << endl;
        c_ld = (t_ld - sum_kahan_ld) - y_ld;    // record the lost low-order bits
        cout << "New low-order bits c_ld = " << fixed << setprecision(32) << c_ld << endl;
        sum_kahan_ld = t_ld;                // update sum
        cout << "Updated sum_kahan_ld = " << fixed << setprecision(32) << sum_kahan_ld << endl;
        cout << "---------------------------------------------" << endl;
    }
    cout << "Final sum_kahan_ld = " << fixed << setprecision(32) << sum_kahan_ld << " sum_ld = " << fixed << setprecision(32) << t_ld << endl;
    cout << string(45, '-') << endl;
    cout << "\n" << endl;
    return 0;
}
```

- Question 2
```C++
#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <random>
#include <iomanip>

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