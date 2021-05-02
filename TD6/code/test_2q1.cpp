#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "ConfusionMatrix.hpp"
#include <cstdlib>
#include <cassert>
using namespace std;

/** @file
 * Test suite for the KNN class.
 * This executable will test the ConfusionMatrix class.
*/

// Prints correct usage
void msgleave(char* argv[]) {
        std::cout << "Usage: " << argv[0] << std::endl;
}

int main(int argc, char* argv[]) {
    // Tests correct usage
	if (argc != 1) {
		msgleave(argv);
		return 1;
	}

    // ConfusionMatrix
    ConfusionMatrix confusion_matrix;
    
    // Test initialization
    assert(("[NOK]", confusion_matrix.GetFN() == 0));
    assert(("[NOK]", confusion_matrix.GetFP() == 0));
    assert(("[NOK]", confusion_matrix.GetTN() == 0));
    assert(("[NOK]", confusion_matrix.GetTP() == 0));

    // Test good entry
    confusion_matrix.AddPrediction(0,0);
    assert(("[NOK]", confusion_matrix.GetFN() == 0));
    assert(("[NOK]", confusion_matrix.GetFP() == 0));
    assert(("[NOK]", confusion_matrix.GetTN() == 1));
    assert(("[NOK]", confusion_matrix.GetTP() == 0));

    // Test good entry
    confusion_matrix.AddPrediction(1,1);
    assert(("[NOK]", confusion_matrix.GetFN() == 0));
    assert(("[NOK]", confusion_matrix.GetFP() == 0));
    assert(("[NOK]", confusion_matrix.GetTN() == 1));
    assert(("[NOK]", confusion_matrix.GetTP() == 1));

    // Test good entry
    confusion_matrix.AddPrediction(0,1);
    assert(("[NOK]", confusion_matrix.GetFN() == 0));
    assert(("[NOK]", confusion_matrix.GetFP() == 1));
    assert(("[NOK]", confusion_matrix.GetTN() == 1));
    assert(("[NOK]", confusion_matrix.GetTP() == 1));

    // Test good entry
    confusion_matrix.AddPrediction(1,0);
    assert(("[NOK]", confusion_matrix.GetFN() == 1));
    assert(("[NOK]", confusion_matrix.GetFP() == 1));
    assert(("[NOK]", confusion_matrix.GetTN() == 1));
    assert(("[NOK]", confusion_matrix.GetTP() == 1));

    // Test metrics
    assert(("[NOK]", confusion_matrix.precision() == 0.5));
    assert(("[NOK]", confusion_matrix.f_score() == 0.5));
    assert(("[NOK]", confusion_matrix.false_alarm_rate() == 0.5));
    assert(("[NOK]", confusion_matrix.detection_rate() == 0.5));
    assert(("[NOK]", confusion_matrix.error_rate() == 0.5));

    // Print
    confusion_matrix.PrintEvaluation();
    
	return 0;
}
