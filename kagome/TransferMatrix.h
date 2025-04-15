#pragma once

#include <vector>
#include <iostream>
#include <Dense>

#include "State.h"

class TransferMatrix {
    Eigen::MatrixXd elements;
    Eigen::VectorXd maxEigenvector;
    double maxEigenvalue;

public:
    TransferMatrix(std::size_t size) : elements(Eigen::MatrixXd::Zero(size, size)) {}

    double& at(std::size_t i, std::size_t j);
    const double& at(std::size_t i, std::size_t j) const;

    const Eigen::VectorXd& getMaxEigenvector() const;
    double getMaxEigenvalue() const;
    const Eigen::MatrixXd& getElements() const;

    std::size_t size() const;

    void power(unsigned);
    double trace();
    void calcMaxEigen();

    friend std::ostream& operator<<(std::ostream&, const TransferMatrix&);
    friend TransferMatrix operator*(const TransferMatrix&, const TransferMatrix&);
};