#include <bitset>
#include <cctype> //std::isdigit
#include <chrono>
#include <cmath>
#include <iostream>
#include <format>
#include <functional> //std::function
#include <unordered_set>
#include <stdexcept>
#include <string>
#include <string_view>

#include <Dense>

#include "Polynomial.h"
#include "State.h"
#include "Timer.h"

template<typename _Scalar>
concept EigenScalar = requires {
    Eigen::Matrix<_Scalar, -1, -1>();
};

double maxEigen(Eigen::MatrixXd mat, unsigned n) {
    Eigen::VectorXd testVect = Eigen::VectorXd::Random(mat.cols());
    testVect /= testVect.norm();
    double eigenvalue = 0;

    for (unsigned i = 0; i < n; i++) {
        testVect = mat * testVect;
        eigenvalue = testVect.norm();

        testVect /= eigenvalue;
    }

    return eigenvalue;
}

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

template<EigenScalar _Scalar>
Eigen::Matrix<_Scalar, -1, -1> genTransfer(std::function<_Scalar(State, State)> countFunc, unsigned dim) {
    Eigen::Matrix<_Scalar, -1, -1> resultMatrix = Eigen::Matrix<_Scalar, -1, -1>::Zero(static_cast<int>(std::pow(2, dim)), static_cast<int>(std::pow(2, dim)));
    std::size_t i = 0, maxj = 0, j;
    State inState(dim), outState(dim);

    for (std::size_t m = 0; m < dim + 1; m++) {
        inState.begin(m);

        while (true) {
            j = maxj;
            outState.begin(m);

            do {
                resultMatrix(i, j) = countFunc(inState, outState);
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

Polynomial countSqPoly(State inState, State outState) {
    int exchangeIndex = -1;
    bool statesIn;
    size_t xCount = 0, yCount = 0;
    
    for (size_t i = 0; i < inState.size() && exchangeIndex == -1; i++) {
        if (inState.at(i) != outState.at(i)) {
            exchangeIndex = i;
            statesIn = inState.at(i);
        }
    }

    if (exchangeIndex == -1) {
        Polynomial result(inState.size());
        result.coefficients(inState.size(), 0) = 1;
        result.coefficients(0, inState.size()) = 1;

        return result;
    }

    for (size_t i = 0; i < exchangeIndex; i++) {
        (inState.at(i) == statesIn ? xCount : yCount)++; //funny syntax
    }

    for (size_t i = exchangeIndex + 1; i < inState.size(); i++) {
        if (inState.at(i) == outState.at(i)) {
            (inState.at(i) == inState.at(exchangeIndex) ? yCount : xCount)++;
        }
        else {
            if (inState.at(i) == inState.at(exchangeIndex)) return Polynomial(0);

            exchangeIndex = i;
        }
    }

    if (inState.at(exchangeIndex) == statesIn) return Polynomial(0);

    Polynomial result(inState.size());
    result.coefficients(xCount, yCount) = 1;
    
    return result;
}

Polynomial countDcPoly(State inState, State outState) {
    if (inState.size() % 2 != 0)
        throw std::invalid_argument("Wymiar musi byc parzysty, aby moc obliczyc macierz dla nieodzialujacych atomow");

    size_t xCount = 0, yCount = 0;

    if (inState.at(0) == outState.at(outState.size() - 1)) {
        if (inState.at(1) == outState.at(0)) {
            if (inState.at(0) == inState.at(1)) {
                xCount++;
            }
            //if we get here its a z-type node
        }
        else {
            return Polynomial(0);
        }
    }
    else if (inState.at(0) == outState.at(0) && inState.at(1) == outState.at(outState.size() - 1)) {
        yCount++;
    }
    else {
        return Polynomial(0);
    }

    for (size_t i = 2; i < inState.size(); i += 2) {
        if (inState.at(i) == outState.at(i - 1)) {
            if (inState.at(i + 1) == outState.at(i)) {
                if (inState.at(i) == inState.at(i + 1)) {
                    xCount++;
                }
                //if we get here its a z-type node
            }
            else {
                return Polynomial(0);
            }
        }
        else if (inState.at(i) == outState.at(i) && inState.at(i + 1) == outState.at(i - 1)) {
            yCount++;
        }
        else {
            return Polynomial(0);
        }
    }

    Polynomial result(inState.size() / 2);
    result.coefficients(xCount, yCount) = 1;

    return result;
}

template<EigenScalar _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols>& power(Eigen::Matrix<_Scalar, _Rows, _Cols>& mat, unsigned n) {
    Eigen::Matrix<_Scalar, _Rows, _Cols> original = mat;

    for (unsigned i = static_cast<unsigned>(std::log2(n)); i > 0; i--) {
        mat *= mat;
        if (n & (1 << i - 1)) {
            mat *= original;
        }
    }

    return mat;
}

void runCalcTransfer(std::string type, unsigned nrows, unsigned ncols) {
    Timer stopwatch;
    unsigned nodeCount;
    double stateCount;

    if (type.substr(type.size() - 4) != "Poly") {
        Eigen::MatrixXd tMat;

        if (type == "square") {
            stopwatch.start();
            tMat = genTransfer(std::function(countSqStates), nrows);
            stopwatch.end();
            nodeCount = nrows * ncols;
        }
        else if (type == "decoupled") {
            stopwatch.start();
            tMat = genTransfer(std::function(countDcStates), nrows);
            stopwatch.end();
            nodeCount = nrows / 2 * ncols;
        }
        else if (type == "kagome") {
            stopwatch.start();
            Eigen::MatrixXd square = genTransfer(std::function(countSqStates), nrows);
            stopwatch.end();
            std::cout << std::format("Square: czas liczenia: {} ms\n", stopwatch.read<Timer::ms>());

            stopwatch.start();
            Eigen::MatrixXd decoupled = genTransfer(std::function(countDcStates), nrows);
            stopwatch.end();
            std::cout << std::format("Decoupled: czas liczenia : {} ms\n", stopwatch.read<Timer::ms>());

            stopwatch.start();
            tMat = square * decoupled;
            stopwatch.end();
            nodeCount = nrows / 2 * 3 * ncols;
        }
        std::cout << std::format("Czas liczenia macierzy transferu: {} ms\n", stopwatch.read<Timer::ms>());

        if (ncols > 1) {
            stopwatch.start();
            power(tMat, ncols);
            stopwatch.end();
            std::cout << std::format("Czas potegowania: {} ms\n", stopwatch.read<Timer::ms>());
        }

        stateCount = tMat.trace();
        std::cout << std::format("Liczba stanow: {}\n", stateCount);
        std::cout << std::format("Entropia resztkowa na czasteczke: {}k_B\n", std::log(stateCount) / nodeCount);
    }
    else {
        Eigen::MatrixXPoly polyMat;

        if (type == "kagomePoly") {
            stopwatch.start();
            Eigen::MatrixXPoly square = genTransfer(std::function(countSqPoly), nrows);
            stopwatch.end();
            std::cout << std::format("Square: czas liczenia: {} ms\n", stopwatch.read<Timer::ms>());

            stopwatch.start();
            Eigen::MatrixXPoly decoupled = genTransfer(std::function(countDcPoly), nrows);
            stopwatch.end();
            std::cout << std::format("Decoupled: czas liczenia: {} ms\n", stopwatch.read<Timer::ms>());

            stopwatch.start();
            polyMat = square * decoupled;
            stopwatch.end();
        }
        std::cout << std::format("Czas liczenia macierzy transferu: {} ms\n", stopwatch.read<Timer::ms>());

        if (ncols > 1) {
            stopwatch.start();
            power(polyMat, ncols);
            stopwatch.end();
            std::cout << std::format("Czas potegowania: {} ms\n", stopwatch.read<Timer::ms>());
        }

        std::cout << "Wielomian: " << polyMat.trace() << '\n';
    }

    std::cout << '\n';
}

void runCalcEigen(std::string type, unsigned nrows, unsigned powerIterCount = 10) {
    Timer stopwatch;
    double eigenvalue;
    unsigned nodeFactor = nrows;
    Eigen::MatrixXd tMat;

    if (type == "square") {
        stopwatch.start();
        tMat = genTransfer(std::function(countSqStates), nrows);
        stopwatch.end();
    }
    else if (type == "decoupled") {
        stopwatch.start();
        tMat = genTransfer(std::function(countDcStates), nrows);
        stopwatch.end();
    }
    else if (type == "kagome") {
        stopwatch.start();
        Eigen::MatrixXd square = genTransfer(std::function(countSqStates), nrows);
        stopwatch.end();
        std::cout << std::format("Square: czas liczenia: {} ms\n", stopwatch.read<Timer::ms>());

        stopwatch.start();
        Eigen::MatrixXd decoupled = genTransfer(std::function(countDcStates), nrows);
        stopwatch.end();
        std::cout << std::format("Decoupled: czas liczenia : {} ms\n", stopwatch.read<Timer::ms>());

        stopwatch.start();
        tMat = square * decoupled;
        stopwatch.end();
        nodeFactor = nodeFactor * 3 / 2;
    }
    std::cout << std::format("Czas liczenia macierzy transferu: {} ms\n", stopwatch.read<Timer::ms>());

    stopwatch.start();
    eigenvalue = maxEigen(tMat, powerIterCount);
    stopwatch.end();
    std::cout << std::format("Czas liczenia wartosci wlasnej: {} ms\n", stopwatch.read<Timer::ms>());

    std::cout << std::format("Maksymalna wartosc wlasna: {}\n", eigenvalue);
    std::cout << std::format("Graniczna entropia resztkowa na czasteczke: {}k_B\n\n", std::log(eigenvalue) / nodeFactor);
}

std::vector<std::string> split(std::string str) {
    std::string word;
    std::vector<std::string> result;
    size_t start = 0, end = str.find(' ');

    while (end != std::string::npos) {
        word = str.substr(start, end - start);

        if (word != "") {
            result.push_back(word);
        }

        start = end + 1;
        end = str.find(' ', start);
    }

    word = str.substr(start, end - start);

    if (word != "") {
        result.push_back(word);
    }

    return result;
}

/*
* Executes given command if syntax is correct. Stores control flags in bitset flags:
* 1st bit - true if syntax was correct
* 2nd bit - true if program is to be exited
*/
void execute(std::vector<std::string> command, std::bitset<2>& flags) {
    if (command.at(0) == "calcTransfer" || command.at(0) == "t") {
        const std::unordered_set<std::string_view> availableTypes = {
                "square",
                "decoupled",
                "kagome", "kagomePoly"
        };

        if (command.size() < 4 || availableTypes.find(command.at(1)) == availableTypes.end()) {
            flags = 0b0;
            return;
        } //should add a check to make sure the next two arguments are unsigneds

        runCalcTransfer(command.at(1), std::stoul(command.at(2)), std::stoul(command.at(3))); //std::stoui doesnt exist :((
    }
    else if (command.at(0) == "calcEntropyEigen" || command.at(0) == "e") {
        const std::unordered_set<std::string_view> availableTypes = {
                "square",
                "decoupled",
                "kagome"
        };

        if (command.size() < 3 || availableTypes.find(command.at(1)) == availableTypes.end()) {
            flags = 0b0;
            return;
        } //should add a check to make sure the next two arguments are unsigneds

        if (command.size() > 3) {
            runCalcEigen(command.at(1), std::stoul(command.at(2)), std::stoul(command.at(3))); 
        }
        else {
            runCalcEigen(command.at(1), std::stoul(command.at(2))); //std::stoui doesnt exist :((
        }
    }
    else if (command.at(0) == "exit") {
        flags = 0b11;
        return;
    }
    else {
        flags = 0b0;
        return;
    }

    flags = 0b1;
}

int main() {
    std::bitset<2> executeFlags = 0;

    while (!executeFlags.test(1)) {
        std::string commandRaw;
        std::vector<std::string> command;

        std::cout << ">";
        std::getline(std::cin, commandRaw);
        command = split(commandRaw);

        execute(command, executeFlags);

        if (!executeFlags.test(0)) {
            std::cout << "Niepoprawna komenda\n\n";
        }
    }

#ifdef KAGOME_TRANSFER_TEST
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
    kagome.power(dim); //matrix exponention can be imporved - binary exp algorithm + only raise to half of the target power, then use trace of product
    clock.end();
    std::cout << std::format("Matrix exponentiation time: {} ms\n", clock.read<Timer::ms>());

    kagome.calcMaxEigen();

    std::cout << std::format("Number of states: {}", kagome.trace()) << '\n';
    std::cout << std::format("Entropy: {}\n", std::log(kagome.trace()) * 2 / 3 / dim / dim);
    std::cout << std::format("Max eigenvalue: {}\n", kagome.getMaxEigenvalue());

    Eigen::IOFormat fmtVect(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", ";");
    Eigen::IOFormat fmtMat(4, 0, ", ", "\n", "[", "]");
    std::cout << "Corresponding eigenvector: " << kagome.getMaxEigenvector().format(fmtVect);
#endif

#ifdef KAGOME_POLYNOMIAL_TEST

    unsigned dim;
    std::cout << "Podaj wymiar pionowy: ";
    std::cin >> dim;

    Eigen::IOFormat fmtMat(4, 0, ", ", "\n", "[", "]");

    Eigen::MatrixXPoly testSq, testDc;

    testSq = genPolyMatrix(countSqPoly, dim);
    testDc = genPolyMatrix(countDcPoly, dim);
    auto test = testSq * testDc;

    std::cout << test.trace() << '\n';

    //power(test, dim);

    //std::cout << test.trace();
#endif
}