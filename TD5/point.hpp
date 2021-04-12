#pragma once

class point
{
	static int d;

public:
	double *coords;
	int label;

	static bool set_dim (int _d);
	static int get_dim ();

	point ();
	~point ();

	void print();
	double dist (point &q);
};
