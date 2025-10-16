#include "utils_gsl.h"
#include <gsl/gsl_vector.h>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <config file>" << endl;
        return 1;
    }
    string config_name = argv[1];
    map<string, string> config = load_config(config_name);
    if (config.empty()) {
        cerr << "Error: no configuration loaded." << endl;
        return 1;
    }
    cout << "Loaded configuration:" << endl;
    for (const auto& [key, value] : config) {
        cout << key << " = " << value << endl;
    }

    FILE* fx = fopen(config["x_file"].c_str(), "r");
    gsl_vector* x = gsl_vector_alloc(100);
    gsl_vector_fread(fx, x);
    fclose(fx);

    FILE* fy = fopen(config["y_file"].c_str(), "r");
    gsl_vector* y = gsl_vector_alloc(100);
    gsl_vector_fread(fy, y);
    fclose(fy);
    if (x->size != y->size) {
        cerr << "Error: input vectors have different sizes." << endl;
        gsl_vector_free(x);
        gsl_vector_free(y);
        return 1;
    }
    double a = stod(config["a"]);

    gsl_vector* z = gsl_vector_alloc(x->size);
    gsl_vector_memcpy(z, y);
    gsl_vector_axpby(a, x, 1.0, z); // z = a * x + 1.0 * y

    for (size_t i = 0; i < z->size; ++i) {
        cout << "z[" << i << "]=" << gsl_vector_get(z, i) << endl;
    }

    FILE* fd = fopen(config["d_file"].c_str(), "w");
    gsl_vector_fwrite(fd, z);
    fclose(fd);

    gsl_vector_free(x);
    gsl_vector_free(y);
    gsl_vector_free(z);
    return 0;
}