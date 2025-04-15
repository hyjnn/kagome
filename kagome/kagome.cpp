#include <chrono>
#include <cmath>
#include <iostream>
#include <format>
#include <functional>
#include <stdexcept>

#include <Dense>

#include "State.h"
#include "TransferMatrix.h"
#include "Timer.h"

double countSqStates(State inState, State outState) {
    double exchangeType = -1;
    double firstExchange = -1;

    for (std::size_t k = 0; k < inState.size(); k++) {
        if (inState.at(k) != outState.at(k)) {
            if (exchangeType == -1) firstExchange = inState.at(k);
            if (exchangeType == inState.at(k)) return 0;
            exchangeType = inState.at(k);
        }
    }

    if (firstExchange == -1) return 2;
    if (firstExchange == exchangeType) return 0;
    return 1;
}

double countDcStates(State inState, State outState) {
    if (inState.size() % 2 != 0)
        throw std::invalid_argument("Wymiar musi byc parzysty, aby moc obliczyc macierz dla nieodzialujacych atomow");

    if (inState.at(0) != outState.at(outState.size() - 1) || inState.at(1) != outState.at(0))
        if (inState.at(0) != outState.at(0) || inState.at(1) != outState.at(outState.size() - 1))
            return 0;
    for (std::size_t k = 2; k < inState.size(); k += 2)
        if (inState.at(k) != outState.at(k - 1) || inState.at(k + 1) != outState.at(k))
            if (inState.at(k) != outState.at(k) || inState.at(k + 1) != outState.at(k - 1))
                return 0;

    return 1;
}

TransferMatrix genTransfer(std::function<double(State, State)> countFunc, int dim) {
    TransferMatrix resultMatrix(std::pow(2, dim));
    std::size_t i = 0, maxj = 0, j;
    State inState(dim), outState(dim);

    for (std::size_t m = 0; m < dim + 1; m++) {
        inState.begin(m);

        while (true) {
            j = maxj;
            outState.begin(m);

            do {
                resultMatrix.at(i, j) = countFunc(inState, outState);
                j += 1;
            } while (outState.next());

            i += 1;
            if (!inState.next()) {
                maxj = j;
                break;
            }
        }
    }

    return resultMatrix;
}

int main() {
    unsigned dim;
    std::cout << "Podaj wymiar pionowy: ";
    std::cin >> dim;
    Timer clock;
    clock.start();
    TransferMatrix square = genTransfer(countSqStates, dim);
    clock.end();
    std::cout << std::format("Square calculation time: {} ms\n", clock.read<Timer::ms>());

    clock.start();
    TransferMatrix decoupled = genTransfer(countDcStates, dim);
    clock.end();
    std::cout << std::format("Decoupled calculation time: {} ms\n", clock.read<Timer::ms>());

    clock.start();
    TransferMatrix kagome = square * decoupled;
    clock.end();
    std::cout << std::format("Matrix multiplication time: {} ms\n", clock.read<Timer::ms>());

    clock.start();
    kagome.power(dim);
    clock.end();
    std::cout << std::format("Matrix exponentiation time: {} ms\n", clock.read<Timer::ms>());

    kagome.calcMaxEigen();

    std::cout << "Number of states: " << kagome.trace() << '\n';
    std::cout << std::format("Entropy: {}\n", std::log(kagome.trace()) * 2 / 3 / dim / dim);
    std::cout << std::format("Max eigenvalue: {}\n", kagome.getMaxEigenvalue());

    Eigen::IOFormat fmtVect(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", ";");
    Eigen::IOFormat fmtMat(4, 0, ", ", "\n", "[", "]");
    std::cout << "Corresponding eigenvector: " << kagome.getMaxEigenvector().format(fmtVect);
}