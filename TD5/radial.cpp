#include <cmath>

#include "point.hpp"
#include "cloud.hpp"
#include "radial.hpp"
#include <math.h>
#include <iostream>

radial::radial(cloud *data_, double bandwidth_): kernel(data_) {
	bandwidth = bandwidth_;
}

double radial::density(point &p) {
	
	int n = data->get_n();	
	
	double s = 0;
	
	for (int k = 0; k < n; k++) {
		double t = pow(p.dist(data->get_point(k)) / bandwidth, 2);
		s += profile(t);
	}
	
	double d = s / (n * volume() * pow(bandwidth, p.get_dim()));
	
	return d;
}
