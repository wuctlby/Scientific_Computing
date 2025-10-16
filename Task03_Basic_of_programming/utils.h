#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <iostream>
#include <map>
#include <algorithm>
#include "H5Cpp.h"

using namespace std;
using namespace H5;

void write_vector(const string& file_name, const vector<double>& v, const string& name, int dim);
map<string, string> load_config(const string& config_name);
vector<double> load_vector(const string& file_name);
inline void trim(string &s);

inline void trim(string &s) {
    s.erase(s.begin(), find_if(s.begin(), s.end(), [](unsigned char ch){ return !isspace(ch); }));
    s.erase(find_if(s.rbegin(), s.rend(), [](unsigned char ch){ return !isspace(ch); }).base(), s.end());
}

map<string, string> load_config(const string& config_name) {
    map<string, string> config; // no sections, just key-value pairs
    ifstream file(config_name);
    if (!file.is_open()) {
        cerr << "Error: cannot open config file " << config_name << endl;
        return config;
    }
    string line;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#' || line[0] == ';')
            continue; // skip comments and empty lines
        size_t eq_pos = line.find('=');
        if (eq_pos == string::npos)
            continue; // skip malformed lines
        string key = line.substr(0, eq_pos);
        trim(key);
        string value = line.substr(eq_pos + 1);
        trim(value);
        config[key] = value;
    }
    file.close();
    return config;
}

void write_vector(const string& file_name, const vector<double>& v, const string& name, int dim) {
    ofstream file(file_name);
    if (!file.is_open()) {
        cerr << "Error: cannot open output file " << file_name << endl;
        return;
    }
    file << "# " << name << " (dimension: " << dim << ")\n";
    for (const auto& val : v)
        file << std::setprecision(6) << val << "\n";

    const string h5_name = file_name.substr(0, file_name.find_last_of('.')) + ".h5";
    cout << "Also writing HDF5 file: " << h5_name << endl;
    H5File h5file(h5_name, H5F_ACC_TRUNC);
    hsize_t dims[1] = { v.size() };
    DataSpace dataspace(1, dims);
    DataSet dataset = h5file.createDataSet(name, PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(v.data(), PredType::NATIVE_DOUBLE);
    h5file.close();
    file.close();
    cout << "Written " << v.size() << " elements to " << file_name << " and " << h5_name << endl;
}

vector<double> load_vector(const string& file_name) {
    vector<double> vec;
    cout << "Loading vector from " << file_name << endl;
    ifstream file(file_name);
    if (!file.is_open()) {
        cerr << "Error: cannot open vector file " << file_name << endl;
        return vec;
    }
    double val;
    string line;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#')
            continue; // skip comments and empty lines
        try {
            val = stod(line);
            vec.push_back(val);
        } catch (const invalid_argument& e) {
            cerr << "Warning: skipping invalid line in " << file_name << ": " << line << endl;
        }
    }
    file.close();
    return vec;
}