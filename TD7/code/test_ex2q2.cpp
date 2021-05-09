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

const double EPS = 0.0001;

void print2(std::string s, double a, double b) {
    cout << s << " (expected, computed) = (" << a << ", " << b << ") " << endl;
}

double rel_error(double a, double b) {
    return fabs(a - b) / fabs(a);
}

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
    std::string fname(argv[2]);
    bool is_Boston = (fname.find("train_boston_housing.csv") != std::string::npos);
    bool is_wine_red = (fname.find("train_winequality-red.csv") != std::string::npos);
    bool is_wine_white = (fname.find("train_winequality-white.csv") != std::string::npos);


    std::map<int, double> M;
    if (is_Boston) {
        M[3] = 35.0232; M[5] = 31.2295; M[10] = 28.4425;
    }
    if (is_wine_red) {
        M[3] = 0.708514; M[5] = 0.647433; M[10] = 0.641149;
    }
    if (is_wine_white) {
        M[3] = 0.963231; M[5] = 0.84421; M[10] = 0.734569;
    }

    if (M[k] == 0) {
        cout << "\n\nWe cannot validate the error for this dataset with this k" << endl;
        return 0;
    }
    if (rel_error(M[k], mse) < EPS) {
        cout << "[OK]" << endl;
    } else {
        cout << "[NOK]" << endl;
        cout << "ERROR: There is an error in the value of MSE : " << endl;
        print2("  MSE ", M[k], mse);
    }

    return 0;
}
