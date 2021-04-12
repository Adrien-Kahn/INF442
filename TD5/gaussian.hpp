#pragma once

#include "cloud.hpp"
#include "radial.hpp"


class gaussian : public radial {
public:
	gaussian(cloud *data_, double bandwidth_);

	double volume();
	double profile(double t);
	void guess_bandwidth();
};
