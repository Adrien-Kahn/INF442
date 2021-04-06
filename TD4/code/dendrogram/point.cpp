#include "point.hpp" // include the header file

#include <iostream> // for cout, endl for printing
#include <cmath>    // for sqrt, fabs
#include <cassert>  // for assert

point::point() {
    coords = new double[d];

    for (int m = 0; m < d; m++)
        coords[m] = 0.0;
}

point::~point() {
    delete[] coords;
}

void point::print() {
    for (int j = 0; j < d; j++)
        std::cout << coords[j] << ' ';

    std::cout << name;
}

double point::dist(point &q) {
    double sqd = 0.0;

    for (int m = 0; m < d; m++)
        sqd += (coords[m] - q.coords[m]) * (coords[m] - q.coords[m]);

    return sqrt(sqd);
}

// the variable d of the class point is static
// hence, it must be initialised here
int point::d = 0; 

bool point::set_dim(int _d) {
    assert (_d > 0);
    if (d > 0)
        return false;

    d = _d;
    return true;
}

int point::get_dim() {
    return d;
}
