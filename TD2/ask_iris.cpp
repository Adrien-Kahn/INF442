
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "retrieval.hpp"

const int dim = 4;  // dimension (hard-coded)
int n = 10000;      // upper bound on number of data points in R^{dim}
int nt = 1000;      // nt query points



point P[10000];
std::string itypes[10000];


int main(int argc, char* argv[]) {
    
    std::ifstream is("./iris2.data");
    assert(is.is_open());    
    for (int k = 0; k < n ; k++) {
        P[k] = new double[dim]; 
        for (int i = 0; i < dim; i++) {
            is >> P[k][i];
        }
        is >> itypes[k];
    }
    std::cout << "Done reading the iris database! ";
    std::cout << "There are " << n << " observations. \n\n";

    while (true) {
        std::cout << "\nEnter your own measurements (4 space-separated real numbers in [0,10]):\n(type CTRL-D to exit): " << std::endl;
        if (std::cin.peek() == EOF) break;
        point query = new double[dim];
        std::cin >> query[0] >> query[1] >> query[2] >> query[3];
        double distb = DBL_MAX;
        int nn;
        nn = linear_scan(query, dim, P, n);
        // bsearch(tree, query, dim, P, distb, nn);
        std::cout << "\n  -> Your iris is of type " << itypes[nn] << std::endl;
    }
}
