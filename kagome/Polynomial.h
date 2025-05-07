#pragma once

#include <iostream>
#include <vector>

#include <Dense>

//going clockwise, starting in the upper left arrow with (kagome) node positioned horizontally:
//x - 1 1 1 1
//y - 1 0 0 1
//z - 1 1 0 0
class Polynomial {
public:

	Eigen::MatrixXd coefficients;

    Polynomial() = default;
	explicit Polynomial(Eigen::MatrixXd coefList) : coefficients(coefList) {}
	explicit Polynomial(size_t n) : coefficients(Eigen::MatrixXd::Zero(n + 1, n + 1)) {} //This constructor takes the polynomial's degree as its argument

	inline size_t degree() const;

    bool friend operator==(const Polynomial&, const Polynomial&);

	Polynomial& operator+=(const Polynomial&);
	Polynomial& operator-=(const Polynomial&);
	Polynomial& operator*=(const Polynomial&);
	Polynomial& operator*=(double);

	friend Polynomial operator+(Polynomial, const Polynomial&);
	friend Polynomial operator-(Polynomial, const Polynomial&);
	friend Polynomial operator*(const Polynomial&, const Polynomial&);
	friend Polynomial operator*(const Polynomial&, double);
	friend Polynomial operator*(double, const Polynomial&);
	friend std::ostream& operator<<(std::ostream&, const Polynomial&);
};

inline size_t Polynomial::degree() const {
	return coefficients.rows() - 1;
}

namespace Eigen {
    typedef Matrix<Polynomial, Eigen::Dynamic, Eigen::Dynamic> MatrixXPoly;

    template<> struct NumTraits<Polynomial> : GenericNumTraits<Polynomial>
    {
        typedef Polynomial Real;
        typedef Polynomial NonInteger;
        typedef Polynomial Nested;

        static inline Real epsilon() { return Polynomial(0); }
        static inline Real dummy_precision() { return Polynomial(0); }
        static inline int digits10() { return 0; }

        enum {
            IsInteger = 0,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = Eigen::HugeCost,
            AddCost = Eigen::HugeCost,
            MulCost = Eigen::HugeCost
        };
    };
}