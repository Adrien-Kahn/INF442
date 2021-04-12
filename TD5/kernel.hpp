#pragma once

#include "point.hpp"
#include "cloud.hpp"

class kernel
{
protected:
	cloud *data;

public:
	kernel(cloud *data_);

	virtual double density(point &p) = 0;
};
