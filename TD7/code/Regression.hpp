#include "Dataset.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>

#ifndef REGRESSION_HPP
#define REGRESSION_HPP

/** 
	The Regression class is an abstract class that will be the basis of the LinearRegression and KnnRegression classes.
*/
class Regression{
protected:
    /**
      The pointer to a dataset.
    */
	Dataset* m_dataset;
    /**
      The column to do regression on.
    */
	int m_col_regr;
public:
    /**
      The constructor sets private attributes dataset (as a pointer) and the column to do regression on (as an int).
    */
	Regression(Dataset* dataset, int col_regr);
    /**
      The Estimate method is virtual: it will depend on the Regression being of class Linear or Knn.
    */
	virtual double Estimate( const Eigen::VectorXd & x ) const = 0;
    /**
      The getter for m_col_regr
    */
	int GetColRegr() const;
    /**
      The getter for m_dataset
    */
	Dataset* GetDataset() const;
};

#endif //REGRESSION_HPP
