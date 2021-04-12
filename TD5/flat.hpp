#pragma once

#include "cloud.hpp"
#include "radial.hpp"


class flat : public radial {
public:
	flat(cloud *data_, double bandwidth_);

	double volume();
	double profile(double t);
};
