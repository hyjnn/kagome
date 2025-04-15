#include <cstdlib>
#include <format>
#include <numeric>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "TransferMatrix.h"

double& TransferMatrix::at(std::size_t i, std::size_t j) {
    if (i > elements.rows() || j > elements.cols()) throw std::invalid_argument("Accessing out of range element");
    return elements(i, j);
}

const double& TransferMatrix::at(std::size_t i, std::size_t j) const {
    if (i > elements.rows() || j > elements.cols()) throw std::invalid_argument("Accessing out of range element");
    return elements(i, j);
}

const Eigen::VectorXd &TransferMatrix::getMaxEigenvector() const {
    return maxEigenvector;
}

double TransferMatrix::getMaxEigenvalue() const {
    return maxEigenvalue;
}

const Eigen::MatrixXd& TransferMatrix::getElements() const {
    return elements;
}

std::size_t TransferMatrix::size() const {
    return elements.rows();
}

void TransferMatrix::power(unsigned n) {
    Eigen::MatrixXd original = elements;
    for (unsigned i = 1; i < n; i++) {
        elements *= original;
    }
}

double TransferMatrix::trace() {
    return elements.trace();
}

std::ostream& operator<<(std::ostream& stream, const TransferMatrix& mat) {
    return stream << mat.elements;
}

TransferMatrix operator*(const TransferMatrix& A, const TransferMatrix& B) {
    TransferMatrix result(A.size());
    result.elements = A.elements * B.elements;
    return result;
}

void TransferMatrix::calcMaxEigen() {
    Eigen::VectorXd eigenvector(size()), prevVector(size());
    std::vector<double> diffNorms(3);
    double eigenvalue = 0;

    for (auto& x : eigenvector) {
        x = std::rand() % 2;
    }
    eigenvector /= eigenvector.norm();

    do {
        for (std::size_t i = 0; i < diffNorms.size(); i++) {
            prevVector = eigenvector;
            eigenvector = elements * eigenvector;
            eigenvalue = eigenvector.norm();
            eigenvector /= eigenvalue;
            diffNorms.at(i) = (prevVector - eigenvector).norm();
        }
    } while (std::accumulate(diffNorms.begin(), diffNorms.end(), 0.) > 0.01);

    maxEigenvalue = eigenvalue;
    maxEigenvector = eigenvector;
}
