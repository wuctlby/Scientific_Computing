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