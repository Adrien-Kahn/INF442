
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "retrieval.hpp"

const bool dmode = 0;
const bool dverbose = false;
const double deps = 0.00001;

bool test_dist(const std::string fname, const bool verbose, const double eps) {
    std::cout << "Testing the function dist()...\t\t\n";
    std::cout << "  Using data from the file " << fname << std::endl;

    int nb_errors = 0;
    std::ifstream in(fname, std::ios_base::in);
    int N;  // number tests
    int dim;
    in >> N;
    in >> dim;
    std::cout << "  " << N << " pairs of points in dimension " << dim << std::endl;

    point p = new double[dim];
    point q = new double[dim];

    double adist = 0.0;
    double cdist = 0.0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < dim; j++) in >> p[j];
        for (int j = 0; j < dim; j++) in >> q[j];
        in >> adist;
        cdist = dist(p, q, dim);

        if (fabs(adist - cdist) > eps) {
            nb_errors++;
        }
        if (verbose) {
            std::cout << std::endl
                      << "(" << i  << ")  Distance of the points " << std::endl;
            print_point(p, dim);
            print_point(q, dim);
            std::cout << " (estimated,  actual) = ("
                      << adist << ", " << cdist << ")" << std::endl;
        }
    }
    delete[] p;
    delete[] q;

    cout << "   #errors = " << nb_errors << " in distance computations" << endl;
    bool success = (nb_errors == 0);
    std::cout << (success ? "[OK]" : "[NOK]") << std::endl;
    return success;
}

void print_help_msg() {
  std::cout << std::endl
            << "USAGE: To test your implementation of the dist function run the program as follows:" << std::endl
            << std::endl
             << " ./test_distance <mode> <verbose> <epsilon>" << std::endl
            << std::endl
            << "All the parameters are optional but you have give them in that order." << std::endl
            << "<mode>      - could be 0 or 1 - ," << std::endl
            << "\t    0 runs the program with input from file ./my_tests/t_dist.dat" << std::endl
            << "\t    1 runs the program with input from the files in ./tests" << std::endl
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


int main(int argc, char *argv[]) {
    // Read in the arguments
    print_help_msg();

    int arg = 1;
    bool mode = (argc > arg) ? std::stoi(argv[arg]) > 0 : false;
    arg++;
    const bool verbose = (argc > arg) ? std::stoi(argv[arg]) > 0 : dverbose;
    arg++;
    const double eps = (argc > arg) ? std::stod(argv[arg]) : deps;

    if (mode) {
        test_dist("./tests/dist_4.dat", verbose, eps); 
        std::cout << std::endl;
        test_dist("./tests/dist_10.dat", verbose, eps);
        std::cout << std::endl;
        test_dist("./tests/dist_50.dat", verbose, eps);
    } else {
        test_dist("./my_tests/t_dist.dat", verbose, eps);
    }
    return 0;
}
