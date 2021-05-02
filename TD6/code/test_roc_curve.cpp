#include <ctime>
#include <cstdlib>
#include <iostream>
#include <cassert>

#include "Dataset.hpp"
#include "KnnClassification.hpp"
#include "ConfusionMatrix.hpp"

using namespace std;

// Prints correct usage
void msgleave(char* argv[]) {
        std::cout << "Usage: " << argv[0] << " <k> <train_file> <test_file> <num_points> [ <column_for_classification> ]" << std::endl;
}

int main(int argc, char* argv[]) {

    // Checks the number of input arguments
    if (argc < 5) {
        msgleave(argv);
        return 1;
    }

    // Number of nearest neighbours for k-NN
    int k = atoi(argv[1]);
    if (k<1) {
        std::cout<<"[NOK] Invalid value of k."<<std::endl;
        msgleave(argv);
        return 1;
    }
    // Puts train and test files in a Dataset object
    Dataset train_dataset(argv[2]);
    Dataset class_dataset(argv[3]);

    // Number of points to generate the ROC curve
    int num_points = atoi(argv[4]);

    // Checks which column is the class label
    int col_class;
    if (argc == 6) {
        col_class = atoi(argv[5]);
    } else {
        //col_class = train_dataset.getDim()-1;
        col_class = 0;
        // std::cout << "No column specified for classification, assuming first column of dataset ("<< col_class <<")."<<std::endl;
    }

    // Prints dimension
    // train_dataset.Show(false);  // only dim and samples

    // Checks if train and test are same format
    assert(train_dataset.getDim() == class_dataset.getDim());

    // Generates the list of thresholds

    // A confusion matrix for each threshold
    ConfusionMatrix* confusion_matrices = new ConfusionMatrix[num_points];

    // Message printing
    // std::cout<< "Computing k-NN classification (k="<< k << ", classification over column "<< col_class << ")..."<<std::endl;

    // A KNN model
    KnnClassification knn_class(k, &train_dataset, col_class);

    // Prediction for each sample and each threshold
    
    // Adds (1,1) to ROC curve
    cout << 1 <<"\t";
    cout << 1 << endl;
    // Print the couples 'false alarm rate' 'detection rate'
    for (int t = 0; t < num_points; t++) {
        cout << confusion_matrices[t].false_alarm_rate() <<"\t";
        cout << confusion_matrices[t].detection_rate() <<endl;
    }
    // Adds (1,1) to ROC curve
    cout << 0 << "\t";
    cout << 0 << endl;

    // Frees the allocated memory
    delete[] confusion_matrices;
    delete[] thresholds;

    return 0;
}
