#include <cassert>
#include <cfloat>
#include <cmath>
#include <algorithm>

#include "cloud.hpp"

cloud::cloud(int _d, int _nmax)
{
	point::set_dim(_d);

	nmax = _nmax;
	n = 0;

	points = new point[nmax];
}

cloud::~cloud()
{
	delete[] points;
}

int cloud::get_n()
{
	return n;
}

point& cloud::get_point(int i)
{
	return points[i];
}

double cloud::min_in_coord(int m)
{
	double minval = DBL_MAX;

	for(int i = 0; i < n; i++)
		if(points[i].coords[m] < minval)
			minval = points[i].coords[m];

	return minval;
}

double cloud::max_in_coord(int m)
{
	double maxval = DBL_MIN;

	for(int i = 0; i < n; i++)
		if(points[i].coords[m] > maxval)
			maxval = points[i].coords[m];

	return maxval;
}

void cloud::add_point(point &p)
{
	assert(n < nmax);

	for(int m = 0; m < point::get_dim(); m++)
	{
		points[n].coords[m] = p.coords[m];
	}

	n++;
}

void cloud::load(std::ifstream &is) {
	assert(is.is_open());

	// point to read into
	point p;
	p.label = 0;

	// while not at end of file
	while(is.peek() != EOF)
	{
		// read new points
		for(int m = 0; m < point::get_dim(); m++)
		{
			is >> p.coords[m];
		}

		add_point(p);

		// consume \n
		is.get();
	}
}

double cloud::standard_deviation()
{
	point mean;
	for(int m = 0; m < point::get_dim(); m++)
	{
		for(int i = 0; i < n; i++)
			mean.coords[m] += points[i].coords[m];
		mean.coords[m] /= n;
	}

	double sum = 0.0;
	for(int i = 0; i < n; i++)
		sum += mean.dist(points[i])*mean.dist(points[i]);

	return std::sqrt(sum/(n-1));
}

