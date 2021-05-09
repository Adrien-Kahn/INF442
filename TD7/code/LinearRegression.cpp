#include "LinearRegression.hpp"
#include "Dataset.hpp"
#include "Regression.hpp"
#include<iostream>
#include<cassert>

#include <Eigen/Dense>


LinearRegression::LinearRegression( Dataset* dataset, int col_regr ) 
: Regression(dataset, col_regr) {
	SetCoefficients();
}

LinearRegression::~LinearRegression() {
	m_beta->resize(0);
	delete m_beta;
}

void LinearRegression::SetCoefficients() {

	int d = m_dataset->GetDim();
	int n = m_dataset->GetNbrSamples();
	
	Eigen::MatrixXd X(n, d);
	Eigen::VectorXd y(n);
	
	for (int i = 0; i < n; i++) {
		
		X(i, 0) = 1;
		y(i) = m_dataset->GetInstance(i)[m_col_regr];
		
		int k = 1;
		for (int j = 0; j < d; j++) {
			
			if (j != m_col_regr) {
				X(i, k) = m_dataset->GetInstance(i)[j];
				k++;
			}
			
		}
	
	}
	
	// We want to solve Xt * X * beta = Xt * y
	
	Eigen::VectorXd Xty = X.transpose() * y;
	Eigen::MatrixXd XtX = X.transpose() * X;
	
	m_beta = new Eigen::VectorXd;
	*m_beta = XtX.ldlt().solve(Xty);

}

const Eigen::VectorXd* LinearRegression::GetCoefficients() const {
	if (!m_beta) {
		std::cout<<"Coefficients have not been allocated."<<std::endl;
		return NULL;
	}
	return m_beta;
}

void LinearRegression::ShowCoefficients() const {
	if (!m_beta) {
		std::cout<<"Coefficients have not been allocated."<<std::endl;
		return;
	}
	
	if (m_beta->size() != m_dataset->GetDim()) {  // ( beta_0 beta_1 ... beta_{d} )
		std::cout<< "Warning, unexpected size of coefficients vector: " << m_beta->size() << std::endl;
	}
	
	std::cout<< "beta = (";
	for (int i=0; i<m_beta->size(); i++) {
		std::cout<< " " << (*m_beta)[i];
	}
	std::cout << " )"<<std::endl;
}

void LinearRegression::SumsOfSquares( Dataset* dataset, double& ess, double& rss, double& tss ) const {
	assert(dataset->GetDim()==m_dataset->GetDim());

	int d_test = dataset->GetDim();
	int n_test = dataset->GetNbrSamples();
	
	Eigen::MatrixXd X(n_test, d_test);
	Eigen::VectorXd y(n_test);
	
	for (int i = 0; i < n_test; i++) {
		
		X(i, 0) = 1;
		y(i) = dataset->GetInstance(i)[m_col_regr];
		
		int k = 1;
		for (int j = 0; j < d_test; j++) {
			
			if (j != m_col_regr) {
				X(i, k) = dataset->GetInstance(i)[j];
				k++;
			}
			
		}
	
	}
	
	Eigen::VectorXd yhat = X * (*m_beta);
	Eigen::VectorXd ybar = Eigen::VectorXd::Constant(n_test, y.mean());
	
	ess = (yhat - ybar).dot(yhat - ybar);
	rss = (yhat - y).dot(yhat - y);
	tss = (y - ybar).dot(y - ybar);
	
}

double LinearRegression::Estimate( const Eigen::VectorXd & x ) const {
	
	int d = m_dataset->GetDim();
	double yhat = (*m_beta)[0];
	
	int k = 0;
	for (int j = 1; j < d; j++) {
	
		if (k == m_col_regr) k++;
		yhat += (*m_beta)[j] * x[k];
		k++;
	
	}
	
	return yhat;

}
