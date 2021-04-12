#pragma once

#include <fstream>

#include "point.hpp"

class cloud {
	int n;
	int nmax;	// maximum possible number of points
	point *points;

public:
	cloud(int _d, int _nmax);
	~cloud();

	// getters
	int get_n();
	point& get_point(int i);

	// helper functions
	void add_point(point &p);
	void load(std::ifstream &is);

	// computations
	double min_in_coord(int m);
	double max_in_coord(int m);
	double standard_deviation();
};
