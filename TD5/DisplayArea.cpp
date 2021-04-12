#include "DisplayArea.hpp"

bool DisplayArea::on_draw(const Cairo::RefPtr<Cairo::Context> &cr)
{
	Gtk::Allocation allocation = get_allocation();
	const int width = allocation.get_width();
	const int height = allocation.get_height();

	const int x_coord = 0;
	const int y_coord = 1;

	// find min and max on each axis
	double x_min = c->min_in_coord(x_coord);
	double x_max = c->max_in_coord(x_coord);
	double y_min = c->min_in_coord(y_coord);
	double y_max = c->max_in_coord(y_coord);

	// plot all points
	for(int i = 0; i < c->get_n(); i++)
	{
		cr->save();
		cr->set_source_rgba(0.0, 0.0, 0.0, 1.0); 
		cr->arc((c->get_point(i).coords[x_coord]-x_min)*width/(x_max-x_min), (c->get_point(i).coords[y_coord]-y_min)*height/(y_max-y_min), 1.0, 0.0, 2.0 * M_PI); // full circle
		cr->fill_preserve();
		cr->restore();
		cr->stroke();
	}

	return true;
}
