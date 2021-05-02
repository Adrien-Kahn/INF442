#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "RandomProjection.hpp"
#include <cstdlib>
using namespace std;

/** @file
 * Test suite for Exercise 3 question 1.
 * This executable will try to generate a d x l projection matrix
*/

// Prints correct usage
void msgleave(char* argv[]) {
        std::cout << "Usage: " << argv[0] << " <d> <l>" << std::endl;
}

int main(int argc, char* argv[]) {
    // Tests correct usage
	if (argc < 2) {
		msgleave(argv);
		return 1;
	}
	
    // Tests value of k (should be > 1)
	int d=atoi(argv[1]);
	if (d<1) { 
		std::cout<<"[NOK] Invalid value of d."<<std::endl;
		msgleave(argv);
		return 1;
	}

    // Tests value of l (should be < d)
	int l=atoi(argv[2]);
	if (l>d) { 
		std::cout<<"[NOK] Invalid value of l."<<std::endl;
		msgleave(argv);
		return 1;
	}

    Eigen::MatrixXd random_gaussian = RandomProjection::RandomGaussianMatrix(d, l);

    std::cout << "Has d rows?" << std::endl;
    assert(("[NOK]", random_gaussian.rows() == d));
    std::cout << "[OK]" << std::endl;
    std::cout << "Has l columns?" << std::endl;
    assert(("[NOK]", random_gaussian.cols() == l));
    std::cout << "[OK]" << std::endl;
    std::cout << "Estimated mean should be close to 0:" << std::endl;
    std::cout << random_gaussian.mean() << std::endl;
    std::cout << "Estimated variance should be close to: " << 1.0/l << std::endl;
    std::cout << 1.0 / (random_gaussian.size() - 1) * random_gaussian.squaredNorm() << std::endl;

    return 0;
}
