
#include "KnnRegression.hpp"
#include<iostream>
#include <ANN/ANN.h>


KnnRegression::KnnRegression( int k, Dataset* dataset, int col_regr)
: Regression(dataset, col_regr) {
	m_k = k;
 	// TODO Exercise 2.1
}

KnnRegression::~KnnRegression() {
	// TODO Exercise 2.1
}

double KnnRegression::Estimate(const Eigen::VectorXd & x) const {
	assert(x.size()==m_dataset->GetDim()-1);
	// TODO Exercise 2.2
}

int KnnRegression::GetK() const {
	return m_k;
}

ANNkd_tree* KnnRegression::GetKdTree() const {
	return m_kdTree;
}
