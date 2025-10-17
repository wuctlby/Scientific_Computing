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
