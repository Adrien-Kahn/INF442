#include <iostream>
#include <cstdlib>

#include <gtkmm/application.h>
#include <gtkmm/window.h>
#include <gtkmm/drawingarea.h>

#include "point.hpp"
#include "cloud.hpp"
#include "knn.hpp"
#include "KernelArea.hpp"

int main(int argc, char **argv)
{
	int k = 10;
	double total_volume = 1.0;

	if(argc != 2 && argc != 3 && argc != 4)
	{
		std::cerr << "Usage: " << argv[0] << " datafile [k] [total_volume]" << std::endl;
		std::exit(1);
	}

	if(argc >= 3)
	{
		k = atoi(argv[2]);
	}

	if(argc >= 4)
	{
		total_volume = atof(argv[3]);
	}

	std::ifstream is(argv[1]);

	cloud c(2, 130000);
	c.load(is);

	knn ker(&c, k, total_volume);
	point p;

	// launch graphical interface
	Glib::RefPtr<Gtk::Application> app = Gtk::Application::create();

	Gtk::Window win;
	win.set_title("Density Estimation with Knn Kernel");
	win.set_default_size(800, 800);

	KernelArea area(&c, &ker);
	win.add(area);
	area.show();

	return app->run(win);
}
