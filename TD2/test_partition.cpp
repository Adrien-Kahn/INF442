
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

bool test_partition(const std::string fname, bool verbose, const double eps) {
    std::cout << "Testing the function partition()...\t\t\n";
    std::cout << "  Using data from the file " << fname << std::endl;

    std::ifstream in(fname, std::ios_base::in);
    int N;    // number of points
    int dim;  // dimension
    int T;    // number of tests
    in >> N;
    in >> dim;
    std::cout << "  " << N << " #points in dimension " << dim << std::endl;
    in >> T;
    std::cout << "  There are " << T << " tests for computing the median and partition \n";

    point P[10000];
    point Q[10000];  // A copy of the P

    int nb_errors_m = 0;  // number of errors in median computations
    int nb_errors_p = 0;  // number of errors in partition computations

    for (int i = 0; i < N; i++) {
        P[i] = new double[dim];
        Q[i] = new double[dim];
        for (int j = 0; j < dim; j++) {
            in >> P[i][j];
            Q[i][j] = P[i][j];
        }
    }

    for (int t = 0; t < T; t++) {
        //for (int k=0; k < dim ; k++) { std::cout << P[i][k] << " ";  } ; std::cout << " \n";
        int s, e, c;
        double amedian, cmedian;  // actual and computed median;
        in >> s >> e >> c >> amedian;
        cmedian = compute_median(Q, s, e, c);
        if (fabs(amedian - cmedian) > eps) {
            ++nb_errors_m;
        }
        int idx_median = partition(Q, s, e, c, dim);
        // Check if to the left of median the elements are <=
        for (int j = s; j <= idx_median; j++) {
            if (Q[j][c] > cmedian) {
                nb_errors_p++;
            }
        }
        // Check if to the right of median the elements are >
        for (int j = idx_median + 1; j < e; j++) {
            if (Q[j][c] < cmedian) {
                nb_errors_p++;
            }
        }

        if (verbose) {
            std::cout << std::endl
                      << "(" << t << ")  For the range ["
                      << s << " .. " << e << ") and coord " << c
                      << " the median is (computed,  actual) = ("
                      << cmedian << ", " << amedian << ")" << std::endl;
        }

        for (int i = s; i < e; i++) {
            for (int j = 0; j < dim; j++) {
                Q[i][j] = P[i][j];
            }
        }
    }
    for (int i = 0; i < N; i++) {
        delete[] P[i];
        delete[] Q[i];
    }

    cout << "   #errors = " << nb_errors_m << " in median computations" << endl;
    cout << "   #errors = " << nb_errors_p << " in partition computations" << endl;
    bool success = ((nb_errors_m == 0) && (nb_errors_p == 0));
    std::cout << (success ? "[OK]" : "[NOK]") << std::endl;
    return success;
}

void print_help_msg() {
    std::cout << std::endl
              << "USAGE: To test your implementation of the partition function run the program as follows:" << std::endl
              << std::endl
              << " ./test_partition <mode> <verbose> <epsilon>" << std::endl
              << std::endl
              << "All the parameters are optional but you have give them in that order." << std::endl
              << "<mode>      - could be 0 or 1 - ," << std::endl
              << "\t    0 runs the program with input from file ./my_tests/t_median.dat" << std::endl
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
        test_partition("./tests/partition-01.dat", verbose, eps);
        test_partition("./tests/partition-02.dat", verbose, eps);
        test_partition("./tests/partition-03.dat", verbose, eps);
        test_partition("./tests/partition-04.dat", verbose, eps);
    } else {
        test_partition("./my_tests/t_partition.dat", verbose, eps);
    }

    return 0;
}
