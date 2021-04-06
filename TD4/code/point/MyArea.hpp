#pragma once // single incl.

#include <gtkmm/application.h>
#include <gtkmm/window.h>
#include <gtkmm/drawingarea.h>

#include "cloud.hpp"

// for graphical interface
class MyArea : public Gtk::DrawingArea {
private:
	cloud *c;

public:
	MyArea(cloud *_c) {
		c = _c;
	}

	virtual ~MyArea() {}

protected:
	//Override default signal handler:
	bool on_draw(const Cairo::RefPtr<Cairo::Context> &cr) override;
};

