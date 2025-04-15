#include "State.h"
#include <iostream>
#include <vector>
#include <stdexcept>
#include <string>

const double& State::at(std::size_t i) const {
    return states.at(i);
}

double& State::at(std::size_t i) {
    return states.at(i);
}

std::size_t State::size() const {
    return states.size();
}

std::ostream& operator<<(std::ostream& stream, const State& state) {
    stream << '[';
    for (std::size_t i = 0; i < state.size() - 1; i++) {
        stream << state.at(i) << ", ";
    }
    stream << state.at(state.size() - 1) << ']';

    return stream;
}

void State::begin(std::size_t n) {
    if (n > states.size()) throw std::invalid_argument("State::begin argument cannot exceed state dimension");

    for (std::size_t i = 0; i < n; i++) this->at(i) = 1;

    for (std::size_t i = n; i < states.size(); i++) this->at(i) = 0;
}

bool State::next() {
    int i = states.size() - 2;
    int last = this->at(i + 1);

    while (i >= 0 && this->at(i) == last) i--;

    if (i < 0) return false;
    if (last == 0) {
        this->at(i) = 0;
        this->at(i + 1) = 1;
        return true;
    }
    int trailingOneCount = states.size() - i - 1;

    for (std::size_t k = 0; k < trailingOneCount; k++)
        this->at(states.size() - 1 - k) = 0;

    while (i >= 0 && this->at(i) == 0) i--;

    if (i < 0) return false;

    this->at(i) = 0;
    this->at(i + 1) = 1;
    for (std::size_t k = 0; k < trailingOneCount; k++)
        this->at(i + k + 2) = 1;

    return true;
}