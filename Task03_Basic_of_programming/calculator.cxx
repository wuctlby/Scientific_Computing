#include "utils.h"

using namespace std;
using namespace H5;

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

    vector<double> x = load_vector(config["x_file"]);
    vector<double> y = load_vector(config["y_file"]);
    if (x.empty() || y.empty() || x.size() != y.size()) {
        cerr << "Error: invalid input vectors." << endl;
        return 1;
    }
    double a = stod(config["a"]);

    vector<double> z(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        z[i] = a * x[i] + y[i];
    }

    write_vector(config["d_file"], z, "z", x.size());

    x.clear();
    y.clear();
    z.clear();

    cout << "Reloading vectors from HDF5 files..." << endl;
    string h5_x = config["x_file"].substr(0, config["x_file"].find_last_of('.')) + ".h5";
    H5File h5file_x(h5_x, H5F_ACC_RDONLY);
    H5::DataSet dataset_x = h5file_x.openDataSet("x");
    DataSpace dataspace_x = dataset_x.getSpace();
    hsize_t x_dims[1];
    dataspace_x.getSimpleExtentDims(x_dims, NULL);
    x.resize(x_dims[0]);
    dataset_x.read(x.data(), PredType::NATIVE_DOUBLE);
    cout << "Read x from HDF5 file:" << endl;
    for (size_t i = 0; i < x.size(); ++i) {
        cout << "x[" << i << "]=" << x[i] << endl;
    }
    h5file_x.close();
    
    string h5_y = config["y_file"].substr(0, config["y_file"].find_last_of('.')) + ".h5";
    H5File h5file_y(h5_y, H5F_ACC_RDONLY);
    H5::DataSet dataset_y = h5file_y.openDataSet("y");
    DataSpace dataspace_y = dataset_y.getSpace();
    hsize_t y_dims[1];
    dataspace_y.getSimpleExtentDims(y_dims, NULL);
    y.resize(y_dims[0]);
    dataset_y.read(y.data(), PredType::NATIVE_DOUBLE);
    cout << "Read y from HDF5 file:" << endl;
    for (size_t i = 0; i < y.size(); ++i) {
        cout << "y[" << i << "]=" << y[i] << endl;
    }
    h5file_y.close();

    z.resize(x.size());
    cout << "Computing z = a * x + y with a = " << a << " ..." << endl;
    for (size_t i = 0; i < x.size(); ++i) {
        cout << "x[" << i << "]=" << x[i] << ", y[" << i << "]=" << y[i] << endl;
        z[i] = a * x[i] + y[i];
    }

    write_vector(config["d_file"], z, "hdf5_z", x.size());

    return 0;
}
