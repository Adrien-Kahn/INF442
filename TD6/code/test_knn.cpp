#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <ANN/ANN.h>
#include "Dataset.hpp"
#include "KnnClassification.hpp"
#include "ConfusionMatrix.hpp"
#include <cstdlib>
#include <cassert>
using namespace std;

/** @file
 * Test suite for the KNN class.
 * This executable will put the two provided CSV files (train and test) in objects of class Dataset, perform kNN classification on the provided column with the user-provided k, and print the resulting test set ConfusionMatrix.
*/

// Prints correct usage
void msgleave(char* argv[]) {
        std::cout << "Usage: " << argv[0] << " <k> <train_file> <test_file> [ <column_for_classification> ]" << std::endl;
}

int main(int argc, char* argv[]) {
    // Tests correct usage
	if (argc < 4) {
		msgleave(argv);
		return 1;
	}
	
    // Tests value of k (should be > 1)
	int k=atoi(argv[1]);
	if (k<1) { 
		std::cout<<"[NOK] Invalid value of k."<<std::endl;
		msgleave(argv);
		return 1;
	}
	
    // Puts train and test files in a Dataset object
	Dataset train_dataset(argv[2]);
	Dataset class_dataset(argv[3]);
    
    // Checks which column is the class label
    int col_class;
    if (argc == 5) {
    	col_class = atoi(argv[4]);
    } else {
		//col_class = train_dataset.getDim()-1;
        col_class = 0;
		std::cout<< "No column specified for classification, assuming first column of dataset ("<< col_class <<")."<<std::endl;
    }

    // Prints dimension
	train_dataset.Show(false);  // only dim and samples

    // Checks if train and test are same format
	assert(train_dataset.getDim() == class_dataset.getDim());
	
    // Classification
    std::cout<< "Computing k-NN classification (k="<< k << ", classification over column "<< col_class << ")..."<<std::endl;
    KnnClassification knn_class(k, &train_dataset, col_class);

    // ConfusionMatrix
    ConfusionMatrix confusion_matrix;
	
    // Starts predicting
 	std::cout<< "Prediction and Confusion Matrix filling" <<std::endl;
    clock_t t = clock();
    for (int i=0; i<class_dataset.getNbrSamples(); i++) {
        std::vector<double> sample = class_dataset.getInstance(i);
        ANNpoint query = annAllocPt(class_dataset.getDim()-1, 0);
        int true_label = -1;  // To not leave it uninitialized + will error in AddPrediction
        for (int j=0, j2=0; j<train_dataset.getDim()-1 && j2<train_dataset.getDim(); j++, j2++) {
            if (j==col_class && j2==col_class) {
                true_label = sample[j2];
                j--;
                continue;
            }
            query[j] = sample[j2];
        }
        int predicted_label = knn_class.Estimate(query);
        annDeallocPt(query);
        confusion_matrix.AddPrediction(true_label, predicted_label);
    }
    
    t = clock() - t;

    cout <<endl
         <<"execution time: "
         <<(t*1000)/CLOCKS_PER_SEC
         <<"ms\n\n";

    confusion_matrix.PrintEvaluation();
    
	return 0;
}
