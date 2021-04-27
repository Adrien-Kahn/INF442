#pragma once

#include "point.hpp"
#include "cloud.hpp"
#include "kernel.hpp"

class knn : public kernel {

public:

	int k;
	double V;

	knn(cloud *data, int k_, double V_);
	double density(point &p);

};
