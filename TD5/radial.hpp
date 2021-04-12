#pragma once

#include "point.hpp"
#include "cloud.hpp"
#include "kernel.hpp"

class radial : public kernel {
	
public:
	
	double bandwidth;
	
	radial(cloud *data_, double bandwidth_);
	double density(point &p);
	
	virtual double volume() = 0;
	virtual double profile(double t) = 0;

};
