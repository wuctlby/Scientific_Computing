#include "utils_gsl.h"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <dimension>;" << endl;
        cerr << "       " << argv[1] << " <output directory> default ./" << endl;
        return 1;
    }
    int dim = stoi(argv[1]);
    string output_dir = (argc >= 3) ? argv[2] : "./";
    gsl_vector* x = gsl_vector_calloc(dim);
    gsl_vector* y = gsl_vector_calloc(dim);
    gsl_vector_set_all(x, 0.1);
    gsl_vector_set_all(y, 7.1);
    string fx_name = output_dir + "vector_" + to_string(dim) + "_x.bin";
    string fy_name = output_dir + "vector_" + to_string(dim) + "_y.bin";
    FILE* fx = fopen(fx_name.c_str(), "w");
    gsl_vector_fwrite(fx, x);
    fclose(fx);
    FILE* fy = fopen(fy_name.c_str(), "w");
    gsl_vector_fwrite(fy, y);
    fclose(fy);
    
    return 0;
}