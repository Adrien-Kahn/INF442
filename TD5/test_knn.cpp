#include <iostream>
#include <cmath>
#include <cassert>

#include "point.hpp"
#include "cloud.hpp"
#include "knn.hpp"


int main()
{
	const double eps = 0.001;

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

	knn ker(&c, 3, 2.0);


	p.coords[0] = 0.0;
	p.coords[1] = 0.0;
	assert( std::abs(ker.density(p) - 0.03) < eps );

	p.coords[0] = 42.0;
	p.coords[1] = -1.0;
	assert( std::abs(ker.density(p) - 0.0037677) < eps );

	p.coords[0] = -3.0;
	p.coords[1] = 1.0;
	assert( std::abs(ker.density(p) - 0.0223607) < eps );

	std::cout << "tests for knn passed" << std::endl;

	return 0;
}
