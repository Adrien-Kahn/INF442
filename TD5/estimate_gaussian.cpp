#include <iostream>

#include <gtkmm/application.h>
#include <gtkmm/window.h>
#include <gtkmm/drawingarea.h>

#include "point.hpp"
#include "cloud.hpp"
#include "gaussian.hpp"
#include "KernelArea.hpp"

int main(int argc, char **argv)
{
	double bandwidth = 1.0;

	if(argc != 2 && argc != 3)
	{
		std::cerr << "Usage: " << argv[0] << " datafile [bandwidth]" << std::endl;
		std::exit(1);
	}

	if(argc == 3)
	{
		bandwidth = atof(argv[2]);
	}

	std::ifstream is(argv[1]);

	cloud c(2, 130000);
	c.load(is);

	gaussian ker(&c, bandwidth);
	point p;

	if(argc < 3)
		ker.guess_bandwidth();

	// launch graphical interface
	Glib::RefPtr<Gtk::Application> app = Gtk::Application::create();

	Gtk::Window win;
	win.set_title("Density Estimation with Gaussian Kernel");
	win.set_default_size(800, 800);

	KernelArea area(&c, &ker);
	win.add(area);
	area.show();

	return app->run(win);
}
