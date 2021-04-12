#include <iostream>
#include <cmath>	// for sqrt, fabs

#include "point.hpp"

using std::cout; using std::endl;

point::point() {
	coords = new double[d];
	for(int m = 0; m < d; m++)
		coords[m] = 0.0;
	label = 0;
}

point::~point()
{
	delete[] coords;
}

void point::print()
{
	cout << coords[0];

	for (int j = 1; j < d; j++)
		cout << ' ' << coords[j];
}

double point::dist(point &q)
{
	double sqd = 0.0;

	for (int m = 0; m < d; m++) 
		sqd += (coords[m]-q.coords[m]) * (coords[m]-q.coords[m]);

	return sqrt(sqd);
}

int point::d = 0;

bool point::set_dim (int _d) {
	if (d > 0)
		return false;

	d = _d;
	return true;
}

int point::get_dim () {
	return d;
}
