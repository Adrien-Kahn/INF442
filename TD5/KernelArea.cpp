#include <cmath>
#include <iostream>

#include "point.hpp"
#include "cloud.hpp"
#include "kernel.hpp"
#include "KernelArea.hpp"


KernelArea::KernelArea(cloud *c_, kernel *k_)
{
	c = c_;
	k = k_;

	amp_factor = 3000.0;
	sq_per_row = 80;

	intensities = new double*[sq_per_row];
	for(int i = 0; i < sq_per_row; i++)
		intensities[i] = new double[sq_per_row];

	// choose which coordinates to use for 2d image
	const int x_coord = 0;
	const int y_coord = 1;

	// find min and max on each axis
	double x_min = c->min_in_coord(x_coord);
	double x_max = c->max_in_coord(x_coord);
	double y_min = c->min_in_coord(y_coord);
	double y_max = c->max_in_coord(y_coord);

	const double sq_area = (x_max-x_min)*(y_max-y_min)/(sq_per_row*sq_per_row);

	// compute squares
	std::cerr << "computing squares... ";
	point sq_center;
	for(int i = 0; i < sq_per_row; i++)
	{
		for(int j = 0; j < sq_per_row; j++)
		{
			if((i+j) % 50 == 0)
			{
				std::cout << '.';
				std::cerr << "\r                                             \rcomputing squares... " << 100*(i*sq_per_row + j)/(sq_per_row*sq_per_row)  << '%';
				std::cerr.flush();
			}

			sq_center.coords[x_coord] =  x_min + (i+0.5)*(x_max-x_min)/sq_per_row;
			sq_center.coords[y_coord] =  y_min + (j+0.5)*(y_max-y_min)/sq_per_row;

			intensities[i][j] = amp_factor*sq_area*k->density(sq_center);
		}
	}
	std::cerr << "\r                                             \rcomputing squares... done" << std::endl;
}

KernelArea::~KernelArea()
{
	for(int i = 0; i < sq_per_row; i++)
		delete[] intensities[i];
	delete[] intensities;

}

bool KernelArea::on_draw(const Cairo::RefPtr<Cairo::Context> &cr)
{
	Gtk::Allocation allocation = get_allocation();
	const int width = allocation.get_width();
	const int height = allocation.get_height();

	// draw squares
	for(int i = 0; i < sq_per_row; i++)
	{
		for(int j = 0; j < sq_per_row; j++)
		{
			cr->save();
			cr->rectangle((double)i*width/sq_per_row, (double)j*height/sq_per_row, (double)width/sq_per_row, (double)height/sq_per_row);
			cr->set_source_rgba(1.0-intensities[i][j], 1.0-intensities[i][j], 1.0-intensities[i][j], 1.0);
			cr->fill_preserve();
			cr->restore();
			cr->stroke();
		}
	}

	return true;
}
