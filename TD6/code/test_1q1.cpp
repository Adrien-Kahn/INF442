#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <ANN/ANN.h>
#include <ANN/ANNperf.h>
#include "Dataset.hpp"
#include "KnnClassification.hpp"
#include <cstdlib>
#include <cassert>
using namespace std;

/** @file
 * Test suite for the KNN class.
 * This executable will put the two provided CSV files (train and test) in objects of class Dataset, 
 * initialize kNN classification on the provided column with the user-provided k.
*/

// Prints correct usage
void msgleave(char* argv[]) {
        std::cout << "Usage: " << argv[0] << " <k> <train_file> [ <column_for_classification> ]" << std::endl;
}

int main(int argc, char* argv[]) {
    // Tests correct usage
	if (argc < 3) {
		msgleave(argv);
		return 1;
	}
	
    // Tests value of k (should be > 1)
	int k=atoi(argv[1]);
	if (k<1) { 
		std::cout<<"Invalid value of k."<<std::endl;
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
		std::cout<< "No column specified for classification, assuming first column of dataset ("<< col_class <<")."<<std::endl;
    }

    // Prints dimension
	train_dataset.Show(false);  // only dim and samples
	
    // Classification
    std::cout<< "Computing k-NN classification (k="<< k << ", classification over column "<< col_class << ")..."<<std::endl;
    KnnClassification knn_class(k, &train_dataset, col_class);

    // Tests
    std::cout << "Correct value of col_class" << std::endl;
    assert(("[NOK]", knn_class.getColClass() == col_class));
    std::cout << "[OK]" << std::endl;
    std::cout << "Correct value of k" << std::endl;
    assert(("[NOK]", knn_class.getK() == k));
    std::cout << "[OK]" << std::endl;
    std::cout << "kdtree" << std::endl;
    ANNkdStats stats;
    knn_class.getKdTree()->getStats(stats);
    std::cout << stats.dim << " : dimension of space (e.g. 1899 for mail_train)" << std::endl;
	std::cout << stats.n_pts << " : no. of points (e.g. 4000 for mail_train)" << std::endl;
	std::cout << stats.bkt_size << " : bucket size" << std::endl;
	std::cout << stats.n_lf << " : no. of leaves (including trivial)" << std::endl;
	std::cout << stats.n_tl << " : no. of trivial leaves (no points)" << std::endl;
	std::cout << stats.n_spl << " : no. of splitting nodes" << std::endl;
	std::cout << stats.n_shr << " : no. of shrinking nodes (for bd-trees)" << std::endl;
	std::cout << stats.depth << " : depth of tree" << std::endl;
	std::cout << stats.sum_ar << " : sum of leaf aspect ratios" << std::endl;
	std::cout << stats.avg_ar << " : average leaf aspect ratio" << std::endl;

	return 0;
}
