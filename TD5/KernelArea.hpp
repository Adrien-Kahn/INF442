#pragma once

#include <gtkmm/application.h>
#include <gtkmm/window.h>
#include <gtkmm/drawingarea.h>

#include "cloud.hpp"
#include "kernel.hpp"

// for graphical interface
class KernelArea : public Gtk::DrawingArea
{
private:
	cloud *c;
	kernel *k;
	int sq_per_row;
	double amp_factor;
	double **intensities;

public:
	KernelArea(cloud *c_, kernel *k_);
	~KernelArea();

protected:
	//Override default signal handler:
	bool on_draw(const Cairo::RefPtr<Cairo::Context> &cr) override;
};

