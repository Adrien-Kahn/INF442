
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "retrieval.hpp"

using std::cout;
using std::cin;
using std::endl;

const bool dmode = 0;
const bool dverbose = 0;
const double deps = 0.00001;




bool test_linear_scan(const std::string fname, bool verbose, const double eps) {
    std::cout << "Testing the function defeatist()...\t\t\n";
    std::cout << "  Using data from the file " << fname << std::endl;

    std::ifstream in(fname, std::ios_base::in);

    int N;    // number of points
    int dim;  // dimension
    int T;    // number of tests
    in >> N;
    in >> dim;
    std::cout << "  " << N << " #points in dimension " << dim << std::endl;
    in >> T;
    std::cout << "  There are " << T << " tests \n";

    point P[10000];
    for (int i = 0; i < N; i++) {
        P[i] = new double[dim];
        for (int j = 0; j < dim; j++) {
            in >> P[i][j];    
        }
    }    
    std::cout << "Building kd-tree..." << std::flush;
    node* tree = build(P, 0, N, 0, dim);
    if (tree == NULL) {
        std::cout << "\n Some functions are not implemented " << std::endl; 
            std::cout << "[NOK]" << std::endl;
            return 0;
    }
    std::cout << " done" << std::endl;

    
    point q = new double[dim];
    int nb_errors = 0;
    for (int t = 0; t < T; t++) {
        for (int j = 0; j < dim; j++) { in >> q[j]; }
        int aidx; 
        double adist;
        in >> aidx >> adist;
        
        double cdist = DBL_MAX;
        int cidx = -1;
        defeatist_search(tree, q, dim, P, cdist, cidx);

        if (aidx != cidx) { 
            nb_errors++;
        }
        if (verbose) {
            std::cout << std::endl
                      << "(" << t << ")  For the query point " << std::endl;
            print_point(q, dim);
            std::cout << "the NN has index : (actual, computed) = ( " 
            << aidx << ", " << cidx  << ")" << endl;
         }
    }

    delete[] q;

    cout << "   #errors = " << nb_errors << " in defeatist search" << endl;
    bool success = (nb_errors == 0);
    std::cout << (success ? "[OK]" : "[NOK]") << std::endl;
    return success;
}

void print_help_msg() {
  std::cout << std::endl
            << "USAGE: To test your implementation of the defeatist search algorithm run the program as follows:" << std::endl
            << std::endl
             << " ./test_defeatist <mode> <verbose> <epsilon>" << std::endl
            << std::endl
            << "All the parameters are optional but you have give them in that order." << std::endl
            << "<mode>      - could be 0 or 1 - ," << std::endl
            << "\t    0 runs the program with input from file ./my_tests/t_ls.dat" << std::endl
            << "\t    1 runs the program with input from the files in ./tests/" << std::endl
            << "\t    default value is " << dmode << std::endl
            << std::endl
            << "<verbose>   - could be 0 or 1 - ," << std::endl
            << "\t    0  prints no info on the data set it computes with" << std::endl
            << "\t    1  prints info on the data set it computes with" << std::endl
            << "\t    default value is " << dverbose << std::endl
            << std::endl
            << "<epsilon>   - a double number, which defines the comparison sensitivity" << std::endl
            << "\t    for doubles, default is " << deps << std::endl
            << std::endl << std::endl;
}

int main(int argc, char* argv[]) {
     // Read in the arguments
    print_help_msg();

    int arg = 1;
    bool mode = (argc > arg) ? std::stoi(argv[arg]) > 0 : false;
    arg++;
    const bool verbose = (argc > arg) ? std::stoi(argv[arg]) > 0 : dverbose;
    arg++;
    const double eps = (argc > arg) ? std::stod(argv[arg]) : deps;

    if (mode) {
        test_linear_scan("./tests/def-1000-100-100.dat", verbose, eps);
        test_linear_scan("./tests/def-1000-300-100.dat", verbose, eps);
    } else {
        test_linear_scan("./my_tests/t_def.dat", verbose, eps);
    }

    return 0;
}
