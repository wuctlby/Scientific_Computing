#include "utils.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <dimension>;" << endl;
        cerr << "       " << argv[1] << " <output directory> default ./" << endl;
        return 1;
    }
    int dim = stoi(argv[1]);
    string output_dir = (argc >= 3) ? argv[2] : "./";
    vector<double> x(dim, 0.1);
    vector<double> y(dim, 7.1);
    string fx_name = output_dir + "vector_" + to_string(dim) + "_x.dat";
    string fy_name = output_dir + "vector_" + to_string(dim) + "_y.dat";
    write_vector(fx_name, x, "x", dim);
    write_vector(fy_name, y, "y", dim);
    return 0;
}