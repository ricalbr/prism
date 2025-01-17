#include "../include/prism/utils.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string>

void write_array(const Eigen::VectorXd &arr, const std::string &filename) {
    std::ofstream file(filename, std::ofstream::out | std::ofstream::trunc);
    if (file.is_open()) {
        for (const auto &e : arr) {
            file << e << "\n";
        }
        file.close();
    } else {
        std::cerr << "Unable to open file\n";
    }
}
