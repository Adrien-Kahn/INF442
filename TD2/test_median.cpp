
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "retrieval.hpp"

const bool dmode = 0;
const bool dverbose = 0;
const double deps = 0.00001;

using std::cout;
using std::endl;

bool test_median(const std::string fname, bool verbose, const double eps) {
    std::cout << "Testing the function compute_median()...\t\t\n";
    std::cout << "  Using data from the file " << fname << std::endl;

    int nb_errors = 0;
    std::ifstream in(fname, std::ios_base::in);
    int N;    // number of points
    int dim;  // dimension
    int T;    // number of tests
    in >> N;
    in >> dim;
    std::cout << "  " << N << " #points in dimension " << dim << std::endl;
    in >> T;
    std::cout << "  There are " << T << " tests for computing the median \n";

    point P[20000];

    for (int i = 0; i < N; i++) {
        P[i] = new double[dim];
        for (int j = 0; j < dim; j++) {
            in >> P[i][j];
        }
        //for (int k=0; k < dim ; k++) { std::cout << P[i][k] << " ";  } ; std::cout << " \n";
    }

    for (int i = 0; i < 2; i++) {
        int s, e, c;
        double amedian, cmedian;  // actual and computed median;
        in >> s >> e >> c >> amedian;
        cmedian = compute_median(P, s, e, c);
        if (fabs(amedian - cmedian) > eps) {
            ++nb_errors;
        }
        if (verbose) {
            std::cout << std::endl
                      << "(" << i + 1 << ")  For the range ["
                      << s << " .. " << e << ") and coord " << c
                      << " the median is (computed,  actual) = ("
                      << cmedian << ", " << amedian << ")" << std::endl;
        }
    }
    for (int i = 0; i < N; i++) {
        delete[] P[i];
    }
    
    cout << "   #errors = " << nb_errors << " in median computations" << endl;
    bool success = (nb_errors == 0);
    std::cout << (success ? "[OK]" : "[NOK]") << std::endl;
    return success;
}

void print_help_msg() {
    std::cout << std::endl
              << "USAGE: To test your implementation of the median function run the program as follows:" << std::endl
              << std::endl
              << " ./test_median <mode> <verbose> <epsilon>" << std::endl
              << std::endl
              << "All the parameters are optional but you have give them in that order." << std::endl
              << "<mode>      - could be 0 or 1 - ," << std::endl
              << "\t    0 runs the program with input from file ./my_tests/t_median.dat" << std::endl
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
              << std::endl
              << std::endl;
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
        test_median("./tests/median-01.dat", verbose, eps);
        test_median("./tests/median-02.dat", verbose, eps);
        test_median("./tests/median-03.dat", verbose, eps);
        test_median("./tests/median-04.dat", verbose, eps);
    } else {
        test_median("./my_tests/t_median.dat", verbose, eps);
    }

    return 0;
}
