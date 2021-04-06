#include <iostream>
#include <cassert>
#include <cfloat>	// for DBL_MAX
#include <fstream>

#include "cloud.hpp"
#include "MyArea.hpp"

using std::ifstream;
using std::cout; using std::endl;

int main(int argc, char **argv) {
    const int d = 4;
    const int nmax = 150;
    const int k = 3;
  
    // construct point cloud
    cloud c(d, nmax, k);
  
    // open data file
    ifstream is("../../data/iris.data");
    assert(is.is_open());
  
    // point to read into
    point p;
  
    // labels to cycle through
    int label = 0;
  
    // while not at end of file
    while (is.peek() != EOF) {
        // read new points
        for (int m = 0; m < d; m++) {
            is >> p.coords[m];
        }
        
        c.add_point(p, label);
  
        label = (label + 1) % k;
  
        // read ground-truth labels
        // unused in normal operation
        std::string next_name;
        is >> next_name;
  		
        // consume \n
        is.get();
    }
  
    // execute k-means algorithm
    cout << "Intracluster variance before k-means: "
         << c.intracluster_variance()
         << endl;
    
    c.kmeans();
  
    cout << "Intracluster variance after k-means: "
         << c.intracluster_variance()
         << endl;
  
  
    // launch graphical interface
    Glib::RefPtr<Gtk::Application> app = Gtk::Application::create(argc, argv, "inf442.td3");
  
    Gtk::Window win;
    win.set_title("TD 3");
    win.set_default_size(400, 400);
  
    MyArea area(&c);
    win.add(area);
    area.show();
  
    return app->run(win);
}
