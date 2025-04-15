#include <algorithm>
#include <stdexcept>
#include <vector>

#include "stateClass.h"

bool operator==(const stateClass& a, const stateClass& b) {
	std::vector<unsigned> compBuffer = b.label;
	
	if (a.dimension != b.dimension) return false;

	for (std::size_t i = 1; i < a.dimension; i++) {
		if (a.label == compBuffer) return true;
		std::rotate(compBuffer.begin(), compBuffer.end() - 1, compBuffer.end());
	}
	
	return false;
}

stateClass::stateClass(unsigned dim, unsigned subdim) {
	initialize(dim, subdim);
}

void stateClass::initialize(unsigned dim, unsigned subdim) {
	if (subdim > dim) throw std::invalid_argument("subdim cannot exceed dim");

	dimension = dim;
	subdimension = subdim;
	label = { subdim, dim - subdim };
}
