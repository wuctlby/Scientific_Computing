#include <fstream>
#include <gsl/gsl_vector.h>
#include <string>
#include <iomanip>
#include <iostream>
#include <map>
#include <algorithm>

using namespace std;

map<string, string> load_config(const string& config_name);
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