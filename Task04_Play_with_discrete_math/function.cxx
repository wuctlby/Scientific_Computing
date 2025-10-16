#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_vector.h>
#include <vector>

using namespace std;

int N{100}; // number of points
double x_inf{0}; // initial value
double x_sup{M_PI_2};

double fun(double x) {
    return gsl_sf_cos(x)*gsl_sf_exp(x);
}

int main (int argc, char *argv[]) {
    if (argv[1]) N = atoi(argv[1]);
    if (argv[2]) x_inf = atoi(argv[2]);
    if (argv[3]) x_sup = atoi(argv[3]);

    double x;
    double range_x = x_sup - x_inf;
    double step_x = range_x / N;
    // double fun = [] (double x) -> double { return gsl_sf_cos(x)*gsl_sf_exp(x); };

    gsl_vector* x_vector = gsl_vector_alloc(N);
    int i=0;
    FILE* dp = fopen("./data_points.dat", "w");
    fprintf(dp, "x\t\t f(x)\n");
    for (x=x_inf; x<=x_sup; x+=step_x) {
        double fx = fun(x);
        gsl_vector_set(x_vector, i, fx);
        ++i;
        fprintf(dp, "%f.4 \t %f.4", x, fx);
    }
    gsl_vector_free(x_vector);
    fclose(dp);
    return 0;
}