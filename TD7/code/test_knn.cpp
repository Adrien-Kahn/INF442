#include <ANN/ANN.h>

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "Dataset.hpp"
#include "KnnRegression.hpp"

using std::cout;
using std::endl;
/** @file
 * Test suite for the KNN class.
 * This executable will put the two provided CSV files (train and test) in objects of class Dataset, perform kNN regression on the provided column with the user-provided k, and print the resulting test set MSE.
*/

void msgleave(char* argv[]) {
    std::cout << "Usage: " << argv[0] << " <k> <train_file> <test_file> [ <column_for_regression> ]" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        msgleave(argv);
        return 1;
    }

    int k = atoi(argv[1]);
    if (k < 1) {
        std::cout << "Invalid value of k." << std::endl;
        msgleave(argv);
        return 1;
    }

    Dataset train_dataset(argv[2]);
    Dataset regr_dataset(argv[3]);

    int col_regr;
    if (argc == 5) {
        col_regr = atoi(argv[4]);
    } else {
        col_regr = train_dataset.GetDim() - 1;
        std::cout << "No column specified for regression, assuming last column of dataset (" << col_regr << ")." << std::endl;
    }

    train_dataset.Show(false);  // only dim and samples

    assert(train_dataset.GetDim() == regr_dataset.GetDim());  // otherwise doesn't make sense

    std::cout << "Computing k-NN regression (k=" << k << ", regression over column " << col_regr << ")..." << std::endl;

    KnnRegression knn_reg(k, &train_dataset, col_regr);

    double mse = 0;
    for (int i = 0; i < regr_dataset.GetNbrSamples(); i++) {
        std::vector<double> sample = regr_dataset.GetInstance(i);
        // extract column for regression
        Eigen::VectorXd query(regr_dataset.GetDim() - 1);
        for (int j = 0, j2 = 0; j < train_dataset.GetDim() - 1 && j2 < train_dataset.GetDim(); j++, j2++) {
            if (j == col_regr && j2 == col_regr) {
                j--;
                continue;
            }
            query[j] = sample[j2];
        }
        double estim = knn_reg.Estimate(query);
        mse += (estim - sample[col_regr]) * (estim - sample[col_regr]) / regr_dataset.GetNbrSamples();
    }

    annClose();

    std::cout << "Mean Square(d) Error (MSE) over test set: " << mse << std::endl;
    

    return 0;
}
