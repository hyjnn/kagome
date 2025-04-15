#pragma once

#include <vector>
#include <iostream>

class State {
    std::vector<double> states;

public:
    State(std::size_t size) : states(std::vector<double>(size)) {}

    const double& at(std::size_t i) const;
    double& at(std::size_t i);
    std::size_t size() const;

    friend std::ostream& operator<<(std::ostream&, const State&);

    void begin(std::size_t n);
    bool next();
};