
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "retrieval.hpp"

int main() {
	
	int const dim = 3;
	int const size = 5;
	
	double a1[] = {7.3, 2.3, 1.9};
	point p1 = a1;
	double a2[] = {0.3, 9.3, 1.1};
	point p2 = a2;
	double a3[] = {3, 4.6, 9.8};
	point p3 = a3;
	double a4[] = {-3.5, 9, 1.4};
	point p4 = a4;
	double a5[] = {-7.8, 9.9, 4};
	point p5 = a5;
	
	point P[] = {p1, p2, p3, p4, p5};

	for (int k = 0; k < size; k++) {
		print_point(P[k], dim);
	}
	
	std::cout << compute_median(P, 1, 5, 0) << std::endl;
	
	return 0;
}
