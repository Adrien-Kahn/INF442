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

    Eigen::VectorXd beta(train_dataset.GetDim());

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

    cout << endl
         << endl;
    double ct_ess, ct_rss, ct_tss;
    tester.SumsOfSquares(&train_dataset, ct_ess, ct_rss, ct_tss);
    std::cout << "Sums of Squares wrt training set:" << std::endl
              << " ESS=" << ct_ess << std::endl
              << " RSS=" << ct_rss << std::endl
              << " TSS=" << ct_tss << std::endl
              << " R2=" << ct_ess / ct_tss << std::endl
              << " MSE=" << ct_rss / train_dataset.GetNbrSamples() << std::endl;
    double cr_ess, cr_rss, cr_tss;
    tester.SumsOfSquares(&regr_dataset, cr_ess, cr_rss, cr_tss);
    std::cout << "Sums of Squares wrt regression set:" << std::endl
              << " ESS=" << cr_ess << std::endl
              << " RSS=" << cr_rss << std::endl
              << " TSS=" << cr_tss << std::endl
              << " MSE=" << cr_rss / regr_dataset.GetNbrSamples() << std::endl;

    return 0;
}
