#pragma once

#include <vector>

class stateClass
{
	std::vector<unsigned> label;
	unsigned dimension, subdimension;

public:
	stateClass(unsigned dim, unsigned subdim);

	void initialize(unsigned dim, unsigned subdim);
	
	friend bool operator==(const stateClass&, const stateClass&);
};

