#include <ANN/ANN.h>
#include <ANN/ANNperf.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Dataset.hpp"
#include "KnnRegression.hpp"
using namespace std;

void print2(std::string s, double a, double b) {
    cout << s << " (expected, computed) = (" << a << ", " << b << ") " << endl;
}

/** @file
 * Test suite for the KNN class.
 * This executable will put the two provided CSV files (train and test) in objects of class Dataset, 
 * initialize kNN regression on the provided column with the user-provided k.
*/

// Prints correct usage
void msgleave(char* argv[]) {
    std::cout << "Usage: " << argv[0] << " <k> <train_file> [ <column_for_regression> ]" << std::endl;
}

int main(int argc, char* argv[]) {
    // Tests correct usage
    if (argc < 3) {
        msgleave(argv);
        return 1;
    }

    // Tests value of k (should be > 1)
    int k = atoi(argv[1]);
    if (k < 1) {
        std::cout << "Invalid value of k." << std::endl;
        msgleave(argv);
        return 1;
    }

    // Puts train and test files in a Dataset object
    Dataset train_dataset(argv[2]);

    // Checks which column is the class label
    int col_class;
    if (argc == 4) {
        col_class = atoi(argv[3]);
    } else {
        //col_class = train_dataset.getDim()-1;
        col_class = 0;
        std::cout << "No column specified for regression, assuming first column of dataset (" << col_class << ")." << std::endl;
    }

    // Prints dimension
    train_dataset.Show(false);  // only dim and samples

    // Regression
    std::cout << "Computing k-NN regression (k=" << k << ", regression over column " << col_class << ")..." << std::endl;
    KnnRegression knn_regr(k, &train_dataset, col_class);

    // Tests

    if (knn_regr.GetColRegr() == col_class) {
        std::cout << "Correct value of col_reg" << std::endl;
        std::cout << "[OK]" << std::endl;
    } else {
        print2("Wrong value of col_reg", knn_regr.GetColRegr(), col_class);
        std::cout << "[NOK]" << std::endl;
    }

    if (knn_regr.GetK() == k) {
        std::cout << "Correct value of k" << std::endl;
        std::cout << "[OK]" << std::endl;
    } else {
        print2("Wrong value of k", knn_regr.GetK(), k);
        std::cout << "[NOK]" << std::endl;
    }

    std::cout << std::endl
              << "Statistics for the ANN kd-tree" << endl;
    ANNkdStats stats;
    knn_regr.GetKdTree()->getStats(stats);
    std::cout << stats.dim << " : dimension of space (e.g. 11 for train_winequality-red)" << std::endl;
    std::cout << stats.n_pts << " : no. of points (e.g. 598 for train_winequality-red)" << std::endl;
    std::cout << stats.bkt_size << " : bucket size" << std::endl;
    std::cout << stats.n_lf << " : no. of leaves (including trivial)" << std::endl;
    std::cout << stats.n_tl << " : no. of trivial leaves (no points)" << std::endl;
    std::cout << stats.n_spl << " : no. of splitting nodes" << std::endl;
    std::cout << stats.n_shr << " : no. of shrinking nodes (for bd-trees)" << std::endl;
    std::cout << stats.depth << " : depth of tree" << std::endl;
    std::cout << stats.sum_ar << " : sum of leaf aspect ratios" << std::endl;
    std::cout << stats.avg_ar << " : average leaf aspect ratio" << std::endl;

    annClose();
    return 0;
}
