#include <iostream>

#include <gtkmm/application.h>
#include <gtkmm/window.h>
#include <gtkmm/drawingarea.h>

#include "point.hpp"
#include "cloud.hpp"
#include "DisplayArea.hpp"

int main(int argc, char **argv)
{
	if(argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " datafile" << std::endl;
		std::exit(1);
	}

	std::ifstream is(argv[1]);

	cloud c(2, 130000);
	c.load(is);

	// launch graphical interface
	Glib::RefPtr<Gtk::Application> app = Gtk::Application::create();

	Gtk::Window win;
	win.set_title("Display Data Set");
	win.set_default_size(800, 800);

	DisplayArea area(&c);
	win.add(area);
	area.show();

	return app->run(win);
}
