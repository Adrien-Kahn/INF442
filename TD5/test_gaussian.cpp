#include <iostream>
#include <cmath>
#include <cassert>

#include "point.hpp"
#include "cloud.hpp"
#include "gaussian.hpp"


int main()
{
	const double eps = 0.01;

	cloud c(7, 5);
	point p;
	p.coords[0] = 1.0;
	p.coords[1] = 2.0;
	p.coords[2] = 3.0;
	p.coords[3] = 4.0;
	p.coords[4] = 5.0;
	p.coords[5] = 6.0;
	p.coords[6] = 7.0;
	c.add_point(p);
	p.coords[0] = 0.0;
	p.coords[1] = -2.0;
	p.coords[2] = 3.0;
	p.coords[3] = 11.0;
	p.coords[4] = -5.0;
	p.coords[5] = 6.0;
	p.coords[6] = 3.0;
	c.add_point(p);

	gaussian ker(&c, 5.00);

	assert( std::abs(ker.volume() - 621.77)/621.77 < eps );

	assert( std::abs(ker.profile(0.5) - 0.778801)/0.778801 < eps );
	assert( std::abs(ker.profile(1.5) - 0.472367)/0.472367 < eps );
	assert( std::abs(ker.profile(2.5) - 0.286505)/0.286505 < eps );

	p.coords[0] = 0;
	p.coords[1] = 0;
	p.coords[2] = 0;
	p.coords[3] = 0;
	p.coords[4] = 0;
	p.coords[5] = 0;
	p.coords[6] = 0;

	assert( std::abs(ker.density(p) - 7.99962e-10)/7.99962e-10 < eps );

	ker.guess_bandwidth();
	assert( std::abs(ker.density(p) - 1.32198e-10)/1.32198e-10 < eps );

	std::cout << "tests for gaussian passed" << std::endl;

	return 0;
}
