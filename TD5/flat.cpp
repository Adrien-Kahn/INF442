#include <cmath>

#include "point.hpp"
#include "flat.hpp"
#include <math.h>

flat::flat(cloud *data_, double bandwidth_): radial(data_, bandwidth) {
}

double flat::volume() {
	double d = (double) data->get_point(0).get_dim();
	return pow(M_PI, d/2) / tgamma((d/2) + 1);
}

double flat::profile(double t) {
	if (t <= 1) return 1;
	return 0;
}
