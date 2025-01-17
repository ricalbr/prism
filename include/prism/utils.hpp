#ifndef PRISM_UTILS_HPP
#define PRISM_UTILS_HPP

#include <Eigen/Dense>
#include <string>

void write_array(const Eigen::VectorXd &arr, const std::string &filename = "c.txt");

#endif // PRISM_UTILS_HPP
