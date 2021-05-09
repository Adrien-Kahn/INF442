#include <math.h>
#include <stdio.h>

#include <cstdlib>
#include <iostream>
#include <string>

#include "Dataset.hpp"
#include "LinearRegression.hpp"

/** @file
 * Test suite for the LinearRegression class.
 * This executable will put the two provided CSV files (train and test) 
 * in objects of class Dataset, perform linear regression on the provided column, 
 * and print the resulting test set MSE.
*/

using std::cout;
using std::endl;

const double EPS = 0.0001;

void print2(std::string s, double a, double b) {
    cout << s << " (expected, computed) = (" << a << ", " << b << ") " << endl;
}

double rel_error(double a, double b) {
    return fabs(a - b) / fabs(a);
}

int main(int argc, const char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <train_file> <test_file> [ <column_for_regression> ]" << std::endl;
        return 1;
    }

    std::cout << std::endl << "Reading training dataset " << argv[1] << " ..." << std::endl << std::endl;

    Dataset train_dataset(argv[1]);
    Dataset regr_dataset(argv[2]);

    int col_regr;
    if (argc == 4) {
        col_regr = atoi(argv[3]);
    } else {
        col_regr = train_dataset.GetDim() - 1;
        std::cout << "No column specified for regression, assuming last column of dataset (" << col_regr << ")." << std::endl;
    }

    std::string fname(argv[1]);

    train_dataset.Show(false);  // only dimensions and samples

    // Checks if train and test are same format
    assert((train_dataset.GetDim() == regr_dataset.GetDim()));  // otherwise doesn't make sense

    cout << endl
         << endl;
    std::cout << "Computing linear regression coefficients (regression over column " << col_regr << ")..." << std::endl;
    LinearRegression tester(&train_dataset, col_regr);

    tester.ShowCoefficients();

    bool is_Boston = (fname.find("train_boston_housing.csv") != std::string::npos);
    bool is_wine_red = (fname.find("train_winequality-red.csv") != std::string::npos);
    bool is_wine_white = (fname.find("train_winequality-white.csv") != std::string::npos);

    Eigen::VectorXd beta(train_dataset.GetDim());
    double tv;                      // true value;
    double ev;                      // estimated value;
 
    if (is_Boston) {
        beta << 30.1835, -0.194652, 0.0440677, 0.0521448, 1.88823, -14.9475, 4.76119,
            0.00262339, -1.30091, 0.46023, -0.0155731, -0.811248, -0.00218155, -0.531514;
        tv = 5;
        ev = 4.24626;
    }
    if (is_wine_red) {
        beta << 15.4472, 0.0335424, -1.01358, -0.277886, 0.0124616, -1.58686,
            0.00670627, -0.0060119, -12.0024, -0.203751, 0.680616, 0.292444;
        tv = 6;
        ev = 5.4228;
    }
    if (is_wine_white) {
        beta << 238.253, 0.141813, -1.90407, -0.0300696, 0.107086, -0.391123,
            0.00603171, -3.67043e-05, -240.576, 1.20104, 0.869698, 0.125772;
        tv = 5;
        ev = 5.88042;
    }

    double dnorm = (beta - *tester.GetCoefficients()).norm();
    double normb = beta.norm();
    if (dnorm / normb <= EPS) {
        cout << "[OK]" << endl;
    } else {
        cout << "[NOK]" << endl;
        cout << "ERROR: There is an error in computing the coefficients: " << endl
             << "  true beta : " << endl
             << beta << endl
             << "  computed beta : " << endl
             << *tester.GetCoefficients() << endl;
    }

    cout << endl
         << endl;
    std::cout << "Testing the Estimate method on the first sample of test_file" << std::endl;
    std::vector<double> first_sample = regr_dataset.GetInstance(0);
    Eigen::VectorXd first_sample_eigen(regr_dataset.GetDim() - 1);

    for (int j = 0, j2 = 0; j < regr_dataset.GetDim() - 1 && j2 < regr_dataset.GetDim(); j++, j2++) {
        if (j == col_regr && j2 == col_regr) {
            j--;
            continue;
        }
        first_sample_eigen(j) = first_sample[j2];
    }

    std::cout << "The true value of y for the first sample of test_file is:" << std::endl;
    std::cout << first_sample[col_regr] << std::endl;
    std::cout << "The estimated value of y for the first sample of test_file is:" << std::endl;
    double cev = tester.Estimate(first_sample_eigen);
    std::cout << cev << std::endl;

    if ((rel_error(tv, first_sample[col_regr]) < EPS) &&
        (rel_error(ev, cev) < EPS)) {
        cout << "[OK]" << endl;
    } else {
        cout << "[NOK]" << endl;
        cout << "ERROR: There is an error in the estimated value: " << endl;
        print2("  true value ", tv, first_sample[col_regr]);
        print2("  estimated value ", ev, cev);
    }

    return 0;
}
