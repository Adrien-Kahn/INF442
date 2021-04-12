#include <iostream>
#include <cmath>
#include <cassert>

#include "point.hpp"
#include "cloud.hpp"
#include "radial.hpp"


class radtest : public radial
{
private:


public:
	radtest(cloud *data_, double bandwidth_) : radial(data_, bandwidth_) {}

	double volume()
	{
		return data->get_n();
	}

	double profile(double t)
	{
		return t;
	}
};


int main()
{
	const double eps = 0.1;

	cloud c(2, 5);
	point p;

	p.coords[0] = 1.0;
	p.coords[1] = 2.0;
	c.add_point(p);
	p.coords[0] = 3.0;
	p.coords[1] = 4.0;
	c.add_point(p);
	p.coords[0] = 0.0;
	p.coords[1] = -1.0;
	c.add_point(p);
	p.coords[0] = 3.0;
	p.coords[1] = 7.0;
	c.add_point(p);
	p.coords[0] = 10.0;
	p.coords[1] = -2.0;
	c.add_point(p);

	radtest ker(&c, 0.42);

	p.coords[0] = 0.0;
	p.coords[1] = 0.0;
	assert( std::abs(ker.density(p) - 248.096) < eps );

	p.coords[0] = 42.0;
	p.coords[1] = -1.0;
	assert( std::abs(ker.density(p) - 9782.45) < eps );

	p.coords[0] = -3.0;
	p.coords[1] = 1.0;
	assert( std::abs(ker.density(p) - 417.779) < eps );

	std::cout << "tests for radial passed" << std::endl;

	return 0;
}
