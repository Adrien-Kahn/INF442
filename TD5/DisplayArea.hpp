#pragma once

#include <gtkmm/application.h>
#include <gtkmm/window.h>
#include <gtkmm/drawingarea.h>

#include "cloud.hpp"

// for graphical interface
class DisplayArea : public Gtk::DrawingArea
{
private:
	cloud *c;

public:
	DisplayArea(cloud *_c)
	{
		c = _c;
	}

	virtual ~DisplayArea() {}

protected:
	//Override default signal handler:
	bool on_draw(const Cairo::RefPtr<Cairo::Context> &cr) override;
};

