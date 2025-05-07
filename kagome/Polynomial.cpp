#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include <Dense>

#include "Polynomial.h"

bool operator==(const Polynomial& p, const Polynomial& q) {
	if (p.degree() != q.degree()) {
		return p.coefficients.isZero() && q.coefficients.isZero();
	}

	return p.coefficients == q.coefficients;
}

Polynomial& Polynomial::operator+=(const Polynomial &p) {
	if (this->degree() != p.degree()) {
		if (*this == Polynomial(0)) {
			this->coefficients = p.coefficients;
		}
		else if (p != Polynomial(0)) {
			throw std::invalid_argument("Cannot add two polynomials of different degrees");
		}
	}
	else {
		this->coefficients += p.coefficients;
	}
	return *this;
}

Polynomial& Polynomial::operator-=(const Polynomial &p) {
	this->coefficients -= p.coefficients;
	return *this;
}

Polynomial operator+(Polynomial p, const Polynomial& q) {
	return p += q;
}

Polynomial operator-(Polynomial p, const Polynomial& q) {
	return p -= q;
}

Polynomial operator*(const Polynomial &p, const Polynomial &q) {
	Polynomial result(p.degree() + q.degree()); //dims are equal to degree + 1

	for (size_t i = 0; i <= result.degree(); i++) {
		for (size_t j = 0; j <= result.degree() - i; j++) { //first two loops select output monomial
			for (size_t k = 0; k <= std::min(i, p.degree()); k++) {
				for (size_t l = 0; l <= std::min(j, p.degree()); l++) { //second two sum over inputs
					if (i - k <= q.degree() && j - l <= q.degree())
						result.coefficients(i, j) += p.coefficients(k, l) * q.coefficients(i - k, j - l); 
				}
			}
		}
	}

	return result;
}

Polynomial operator*(const Polynomial& p, double a) {
	return Polynomial(a * p.coefficients);
}

Polynomial operator*(double a, const Polynomial& p) {
	return Polynomial(a * p.coefficients);
}

Polynomial& Polynomial::operator*=(const Polynomial &p) {
	*this = *this * p;
	return *this;
}

Polynomial& Polynomial::operator*=(double x) {
	this->coefficients *= x;
	return *this;
}

//Matrix form printing
//std::ostream& operator<<(std::ostream& stream, const Polynomial& p) {
//	stream << p.coefficients.format(Eigen::IOFormat(4, 0, ", ", "\n", "[", "]"));
//	return stream;
//}

static std::string monomial(size_t i, size_t j, size_t k) {
	std::string result = "";

	if (i > 1) {
		result += "x^" + std::to_string(i);
	}
	else if (i == 1) result += "x";

	if (j > 1) {
		result += "y^" + std::to_string(j);
	}
	else if (j == 1) result += "y";

	if (k > 1) {
		result += "z^" + std::to_string(k);
	}
	else if (k == 1) result += "z";

	return result;
}

//Polynomial form printing
std::ostream& operator<<(std::ostream& stream, const Polynomial& p) {
	bool isFirst = true;
	std::string textForm = "";

	for (size_t i = 0; i <= p.degree(); i++) {
		for (size_t j = 0; j <= p.degree() - i; j++) {
			if (p.coefficients(i, j) != 0.) {
				if (!isFirst) {
					textForm += " + ";
				}
				else {
					isFirst = false;
				}
				if (p.coefficients(i, j) != 1)
					textForm += std::to_string(p.coefficients(i, j));
				textForm += monomial(i, j, p.degree() - i - j);
			}
		}
	}

	if (isFirst) {
		stream << 0;
	}
	else {
		stream << textForm;
	}

	return stream;
}