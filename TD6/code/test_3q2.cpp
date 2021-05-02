#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "Dataset.hpp"
#include "RandomProjection.hpp"
#include <cstdlib>
using namespace std;

/** @file
*/

// Prints correct usage
void msgleave(char* argv[]) {
        std::cout << "Usage: " << argv[0] << " <k> <projection_dim> <train_file> <sampling> [ <column_for_classification> ]" << std::endl;
}

int main(int argc, char* argv[]) {
    // Tests correct usage
	if (argc < 5) {
		msgleave(argv);
		return 1;
	}
		
    // Puts train file in a Dataset object
	Dataset train_dataset(argv[3]);

    // Tests value of projection_dim (should be > 1 < dimension of dataset)
	int projection_dim=atoi(argv[2]);
	if ((projection_dim<1) | (projection_dim>=train_dataset.getDim()-1)) { 
		std::cout<<"[NOK] Invalid value of projection_dim."<<std::endl;
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

    // Tests value of sampling (should be "Gaussian" or "Rademacher")
	std::string sampling=argv[4];
	if ((sampling!="Gaussian") & (sampling!="Rademacher")) { 
		std::cout<<"[NOK] Invalid value of sampling."<<std::endl;
		msgleave(argv);
		return 1;
	}
    
    // Checks which column is the class label
    int col_class;
    if (argc == 6) {
    	col_class = atoi(argv[5]);
    } else {
		//col_class = train_dataset.getDim()-1;
        col_class = 0;
		std::cout << "No column specified for classification, assuming first column of dataset ("<< col_class <<")."<<std::endl;
    }

    // Prints dimension
	train_dataset.Show(false);  // only dim and samples

	// Random projection
    std::cout << "Performing Random Projection" << std::endl;
    clock_t t_random_projection = clock();
	RandomProjection projection(train_dataset.getDim()-1, col_class, projection_dim, sampling);
    t_random_projection = clock() - t_random_projection;
    std::cout << endl
         <<"Execution time: "
         <<(t_random_projection*1000)/CLOCKS_PER_SEC
         <<"ms\n\n";

	// Trivial tests
	std::cout << "Original dimension is OK?" << std::endl;
	assert(("[NOK]", projection.getOriginalDimension() == train_dataset.getDim() - 1));
    std::cout << "[OK]" << std::endl;
	std::cout << "Classification column is OK?" << std::endl;
	assert(("[NOK]", projection.getColClass() == col_class));
    std::cout << "[OK]" << std::endl;
	std::cout << "Projection dimension is OK?" << std::endl;
	assert(("[NOK]", projection.getProjectionDim() == projection_dim));
    std::cout << "[OK]" << std::endl;
	std::cout << "Type of sampling is OK?" << std::endl;
	assert(("[NOK]", projection.getTypeSample() == sampling));
    std::cout << "[OK]" << std::endl;

	// Tests for projection matrix
	std::cout << "Has d rows?" << std::endl;
    assert(("[NOK]", projection.getProjection().rows() == train_dataset.getDim()- 1));
    std::cout << "[OK]" << std::endl;
    std::cout << "Has l columns?" << std::endl;
    assert(("[NOK]", projection.getProjection().cols() == projection_dim));
    std::cout << "[OK]" << std::endl;
    std::cout << "Estimated mean should be close to 0:" << std::endl;
    std::cout << projection.getProjection().mean() << std::endl;
    std::cout << "Estimated variance should be close to: " << 1.0/projection_dim << std::endl;
    std::cout << 1.0 / (projection.getProjection().size() - 1) * projection.getProjection().squaredNorm() << std::endl;

	// Uncomment this to try your ProjectionQuality bonus function!
    projection.ProjectionQuality(&train_dataset);

	return 0;
}
