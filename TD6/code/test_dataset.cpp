#include "Dataset.hpp"
#include <iostream>
#include <cstdlib>
#include <stdio.h>

/** @file
 * Test suite for the Dataset class.
 *  This executable will put the provided CSV file in an object of class Dataset and invoke its Show method with verbose=true.
 */
int main(int argc, const char * argv[]){
	if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <dataset_file>" << std::endl;
        return 1;
    }

	Dataset dataset(argv[1]);
    
	dataset.Show(false);
	
	return 0;
}
