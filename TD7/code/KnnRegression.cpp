
#include "KnnRegression.hpp"
#include<iostream>
#include <ANN/ANN.h>


KnnRegression::KnnRegression( int k, Dataset* dataset, int col_regr)
: Regression(dataset, col_regr) {
	m_k = k;
	
	int n = dataset->GetNbrSamples();
	int d = dataset->GetDim();
	
	// Builds m_dataPts
	m_dataPts = annAllocPts(n, d - 1);
	
	for (int i = 0; i < n; i++) {
		
		int k = 0;
		for (int j = 0; j < d; j++) {
			
			if (j != col_regr) {
				m_dataPts[i][k] = dataset->GetInstance(i)[j];
				k++;
			}
			
		}
	
	}

	m_kdTree = new ANNkd_tree(m_dataPts, n, d - 1);
	
}

KnnRegression::~KnnRegression() {
	annDeallocPts(m_dataPts);
	delete m_kdTree;
}

double KnnRegression::Estimate(const Eigen::VectorXd & x) const {
	assert(x.size()==m_dataset->GetDim()-1);
	
	int d = m_dataset->GetDim();
	ANNpoint queryPt = annAllocPt(d - 1);
	for (int i = 0; i < d - 1; i++) queryPt[i] = x[i];
	
	ANNidxArray nnIdx = new ANNidx[m_k];
	ANNdistArray dists = new ANNdist[m_k];
	
	m_kdTree->annkSearch(queryPt, m_k, nnIdx, dists);
	
	double yhat = 0;
	
	for (int i = 0; i < m_k; i++) {
		yhat += m_dataset->GetInstance(nnIdx[i])[m_col_regr];
	}
	
	return yhat / ((double) m_k);
	
}

int KnnRegression::GetK() const {
	return m_k;
}

ANNkd_tree* KnnRegression::GetKdTree() const {
	return m_kdTree;
}
