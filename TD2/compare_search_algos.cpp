
#include <cmath>
#include <cstdlib>  // for rand, srand
#include <fstream>
#include <iostream>
#include <string>

using std::cout;
using std::endl;

#include "retrieval.hpp"

const int dim = 4;  // dimension (hard-coded)
int n = 10000;      // upper bound on number of data points in R^{dim}
int nt = 1000;      // nt query points

const int dT = 2;
const double eps = 0.00001;
point P[10000];
std::string itypes[10000];

void print_help_msg() {
    std::cout << std::endl
              << "USAGE: To benchmark your implementation of the search algorithms run the program as follows:" << std::endl
              << std::endl
              << " ./compare_search_algos <tests>" << std::endl
              << std::endl
              << "The parameters is optional ." << std::endl
              << "<test>      - a positive integer- ," << std::endl
              << "\t    the number of tests to perform "<< std::endl
              << "\t    default value is " << dT << std::endl
              << std::endl
              << std::endl;
}

int main(int argc, char* argv[]) {
    // Read in the arguments
    print_help_msg();
    int arg = 1;
    int T = (argc > arg) ? std::stoi(argv[arg]) : dT;

    std::ifstream is("./iris2.data");
    assert(is.is_open());
    for (int k = 0; k < n; k++) {
        P[k] = new double[dim];
        for (int i = 0; i < dim; i++) {
            is >> P[k][i];
        }
        is >> itypes[k];
    }
    std::cout << "Done reading the iris database! ";
    std::cout << "There are " << n << " observations. \n\n";

    std::cout << "Building kd-tree..." << std::flush;
    node* tree = build(P, 0, n, 0, dim);
    if (tree == NULL) {
        std::cout << "\n Some functions are not implemented " << std::endl;
        std::cout << "[NOK]" << std::endl;
        return 0;
    }
    std::cout << " done" << std::endl;

    cout << "We perform " << T << " tests" << endl << endl;
    // Random query points
    point* q = new point[T];
    for (int i = 0; i < T; i++) {
        q[i] = new double[dim];
        for (int j = 0; j < dim; j++) {
            q[i][j] = 10 * (double)rand() / RAND_MAX;
        }
    }

    std::clock_t start, end;

    // Time the linear_scan algorithm
    int idx_ls[T];
    std::cout << "Benchmarking linear scan... " << std::flush;
    start = std::clock();
    for (int i = 0; i < T; i++) {
        idx_ls[i] = linear_scan(q[i], dim, P, n);
    }
    end = std::clock();
    cout << "\n   Total time: " << (((float)(end - start)) )
         << "\t\tavg time: " << (float((end - start) / T))  
         << " us" << endl
         << endl
         << endl;

    // Time the defeatist algorithm
    int idx_def[T];
    int def_errors = 0;
    std::cout << "Benchmarking defeatist..." << std::flush;
    start = std::clock();
    for (int i = 0; i < T; i++) {
        double dist_qP = DBL_MAX;
        defeatist_search(tree, q[i], dim, P, dist_qP, idx_def[i]);
        if (idx_def[i] != idx_ls[i]) {
            def_errors++;
        }
    }
    end = std::clock();
    cout << "\n   Total time: " << (((float)(end - start))  )
         << "\t\tavg time: " << (float((end - start) / T)) 
         << " us" << endl;
    cout << "   #Errors = " << def_errors
         << "\t accuracy: " << 100.0 * (T - def_errors) / T << "% \n";
    cout << endl
         << endl;

    // Time the backtracking algorithm
    int idx_bac[T];
    int bac_errors = 0;
    std::cout << "Benchmarking backtracking..." << std::flush;
    start = std::clock();
    for (int i = 0; i < T; i++) {
        double dist_qP = DBL_MAX;
        backtracking_search(tree, q[i], dim, P, dist_qP, idx_bac[i]);
        if (idx_bac[i] != idx_ls[i]) {
            bac_errors++;
        }
    }
    end = std::clock();
    cout << "\n   Total time: " << (((float)(end - start)))
         << "\t\tavg time: " << (float((end - start) / T)) 
         << " us" << std::endl;
    cout << "   #Errors = " << bac_errors
         << "\t accuracy: " << 100. * (T - bac_errors) / T << "% \n";

    return 0;
}
