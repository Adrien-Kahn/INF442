
#include "KnnClassification.hpp"
#include <iostream>
#include <ANN/ANN.h>

 
KnnClassification::KnnClassification(int k, Dataset* dataset, int col_class)
: Classification(dataset, col_class) {

	m_k = k;
	
	int n = dataset->getNbrSamples();
	int d = dataset->getDim();
	
	// Builds m_dataPts
	m_dataPts = annAllocPts(n, d - 1);
	
	for (int i = 0; i < n; i++) {
		
		int k = 0;
		for (int j = 0; j < d; j++) {
			
			if (j != col_class) {
				m_dataPts[i][k] = dataset->getInstance(i)[j];
				k++;
			}
			
		}
	
	}

	m_kdTree = new ANNkd_tree(m_dataPts, n, d - 1);

}

KnnClassification::~KnnClassification() {
	
	annDeallocPts(m_dataPts);
	delete m_kdTree;
	
}

int KnnClassification::Estimate(const ANNpoint & x, double threshold) {

	ANNidxArray nnIdx = new ANNidx[m_k];
	ANNdistArray dists = new ANNdist[m_k];
	
	m_kdTree->annkSearch(x, m_k, nnIdx, dists);
	
	int count = 0;
	
	for (int i = 0; i < m_k; i++) {
		count += m_dataset->getInstance(nnIdx[i])[m_col_class];
	}
	
	if ( (double) count / (double) m_k > threshold) {
		return 1;
	} else {
		return 0;
	}

}

int KnnClassification::getK() {
    return m_k;
}

ANNkd_tree* KnnClassification::getKdTree() {
    return m_kdTree;
}
