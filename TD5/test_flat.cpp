#include <iostream>
#include <cmath>
#include <cassert>

#include "point.hpp"
#include "cloud.hpp"
#include "flat.hpp"


int main()
{
	const double eps = 0.1;

	cloud c(7, 5);

	flat ker(&c, 5.00);

	assert( std::abs(ker.volume() - 4.72477) < eps );

	assert( std::abs(ker.profile(0.5) - 1.0) < eps );
	assert( std::abs(ker.profile(1.5) - 0.0) < eps );
	assert( std::abs(ker.profile(2.5) - 0.0) < eps );

	std::cout << "tests for flat passed" << std::endl;

	return 0;
}
