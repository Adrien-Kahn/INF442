#include <cmath>
#include <iostream>

#include "point.hpp"
#include "cloud.hpp"
#include "gaussian.hpp"
#include <math.h>

gaussian::gaussian(cloud *data_, double bandwidth_): radial(data_, bandwidth_) {
}

double gaussian::volume() {
	double d = (double) data->get_point(0).get_dim();
	return pow(2 * M_PI, d/2);
}

double gaussian::profile(double t) {
	return exp(-t/2);
}

void gaussian::guess_bandwidth() {
	double n = (double) data->get_n();
	double sigma = data->standard_deviation();
	bandwidth = 1.06 * sigma / pow(n, 0.2);
}
