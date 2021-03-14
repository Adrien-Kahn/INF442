#include <iostream>
#include <cmath>
// Need that for memset
#include <cstring>

/********************
 * Helper functions *
 ********************/

// Read the data matrix from the standard input
void readMatrix (double **matrix, int rows, int columns)
{
    for (int i = 0; i < rows; i++)
    {
        matrix[i] = new double[columns];
        for (int j = 0; j < columns; j++)
            std::cin >> matrix[i][j];
    }
}

// Print an array on the standard output
void printArray (double values[], int length)
{
    for (int i = 0; i < length - 1; i++)
        std::cout << values[i] << " ";

    std::cout << values[length - 1];
}

// Prepare an empty matrix
double **prepareMatrix(int rows, int columns) {
    double **matrix = new double* [rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[columns];
        for (int j = 0; j < columns; j++)
            matrix[i][j] = 0;
    }

    return matrix;
}

/*****************************************************
 * Exercise 3: Mean, variance and standard deviation *
 *****************************************************/

/** 
 * This function computes the mean of the given array of values 
 * 
 * @param values the array with the values
 * @param length the length of the array 'values'
 * @return the mean of the values in the array
 */
double computeMean (double values[], int length)
{
    double sum = 0;

    for (int k = 0; k < length; k++) {
    	sum += values[k];
    }

    return sum / length;
}

/** 
 * This function computes the variance of the given array of values 
 * 
 * @param values the array with the values
 * @param length the length of the array 'values'
 * @return the variance of the values in the array
 */
double computeVariance (double values[], int length)
{
    double mean = computeMean(values, length);
    
    double s = 0;
    
    for (int k = 0; k < length; k++) {
    	s += (values[k] - mean)*(values[k] - mean);
    }
    
    return s/length;
}


/** 
 * This function computes the unbiased sample variance of the given
 * array of values 
 * 
 * @param values the array with the values
 * @param length the length of the array values
 * @return the variance of the values in the array
 */
double computeSampleVariance (double values[], int length)
{
    return ( (double) length / (double) (length - 1)) * computeVariance(values, length);
}

/** 
 * This function computes the standard deviation of the given
 * array of values 
 * 
 * @param values the array with the values
 * @param length the length of the array values
 * @return the variance of the values in the array
 */
double computeStandardDeviation (double values[], int length)
{
    return sqrt(computeVariance(values, length));
}

/** 
 * This function computes the unbiased sample standard deviation
 * of the given array of values 
 * 
 * @param values the array with the values
 * @param length the length of the array values
 * @return the variance of the values in the array
 */
double computeSampleStandardDeviation (double values[], int length)
{
    return sqrt(computeSampleVariance(values, length));
}

/************* Additional tests **************/

void runStudentTests3 (double values[], int length) 
{
    std::cout << "--------------------------------------------------------------------------------" << std::endl
              << "Your own tests could print something here. See the runStudentTests3() function" << std::endl;

    // You can insert any code here — it will be called in the demo mode of Exercise 3

    std::cout << "--------------------------------------------------------------------------------" << std::endl;
}

/*************************************
 * Exercise 4: Working with matrices *
 *************************************/

/**
 * This function prints a rectangular matrix on the standard output, 
 * placing each row on a separate line.  
 * 
 * @param matrix the matrix to print
 * @param rows the number of rows in the matrix
 * @param columns the number of columns
 */
void printMatrix (double **matrix, int rows, int columns) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

/** 
 * This function extracts one row from a data matrix
 * 
 * @param matrix the matrix with the data
 * @param columns the number of columns in the matrix
 * @param index the index of the row to extract
 * @param row the array where the extracted values are to be placed
 */
void getRow (double **matrix, int columns, int index, double row[])
{
	// Clean out completely the row 
	std::memset(row, 0, columns*sizeof(double));
	
	for (int j = 0; j < columns; j++) {
		row[j] = matrix[index][j];
	}
}

/** 
 * This function extracts one column from a data matrix
 * 
 * @param matrix the matrix with the data
 * @param rows the number of rows in the matrix
 * @param index the index of the column to extract
 * @param column the array where the extracted values are to be placed
 */
void getColumn (double **matrix, int rows, int index, double column[])
{
	// Clean out completely the column 
	std::memset(column, 0, rows*sizeof(double));

	for (int i = 0; i < rows; i++) {
		column[i] = matrix[i][index];
	}
}

/************* Additional tests **************/

void runStudentTests4 (double **matrix, int rows, int columns) 
{
    std::cout << "--------------------------------------------------------------------------------" << std::endl
              << "Your own tests could print something here. See the runStudentTests4() function" << std::endl;

    //printMatrix(matrix, rows, columns);

    std::cout << "--------------------------------------------------------------------------------" << std::endl;
}

/****************************************************
 * Exercise 5: Covariance and correlation of arrays *
 ****************************************************/

/**
 * This function computes the covariance of two vectors of data of the same length
 * @param values1 the first vector
 * @param values2 the second vector
 * @param length the length of the two vectors
 * @return the covariance of the two vectors
 */
double computeCovariance(double values1[], double values2[], int length)
{
    double m1 = computeMean(values1, length);
    double m2 = computeMean(values2, length);
    
    double s = 0;
    
    for (int k = 0; k < length; k++) {
    	s += (values1[k] - m1)*(values2[k] - m2);
    }
    
    return s/length;
}

/**
 * This function computes the correlation of two vectors of data of the same length
 * 
 * @param values1 the first vector
 * @param values2 the second vector
 * @param length the length of the two vectors
 * @return the correlation of the two vectors
 */
double computeCorrelation(double values1[], double values2[], int length)
{
	return computeCovariance(values1, values2, length) / (computeStandardDeviation(values1, length) * computeStandardDeviation(values2, length));
}

/************* Additional tests **************/

void runStudentTests5 (double values1[], double values2[], int length) 
{
    std::cout << "--------------------------------------------------------------------------------" << std::endl
              << "Your own tests could print something here. See the runStudentTests5() function" << std::endl;

    // You can insert any code here — it will be called in the demo mode of Exercise 5

    std::cout << "--------------------------------------------------------------------------------" << std::endl;
}

/***************************************************
 * Exercise 6: Covariance and correlation matrices *
 ***************************************************/

/**
 * This function computes the covariance matrix of the matrix provided as argument
 * 
 * @param dataMatrix the input matrix 
 * @param rows the number of rows in the matrix
 * @param columns the number of columns in the matrix
 * @return the covariance matrix
 */
double **computeCovarianceMatrix (double **dataMatrix, int rows, int columns)
{
	// Initialise a square matrix
	double **matrix = prepareMatrix(columns, columns);
	// Prepare temporary storage for columns
	double column1[rows], column2[rows];

	for (int i = 0; i < columns; i++) {
		for (int j = 0; j < columns; j++) {
			getColumn(dataMatrix, rows, i, column1);
			getColumn(dataMatrix, rows, j, column2);
			matrix[i][j] = computeCovariance(column1, column2, rows);
		}
	}
	
	return matrix;
}

/**
 * This function computes the correlation matrix of the matrix provided as argument
 * 
 * @param dataMatrix the input matrix 
 * @param rows the number of rows in the matrix
 * @param columns the number of columns in the matrix
 * @return the correlation matrix
 */
double **computeCorrelationMatrix (double **dataMatrix, int rows, int columns)
{
	// Initialise a square matrix
	double **matrix = prepareMatrix(columns, columns);
	// Prepare temporary storage for columns
	double column1[rows], column2[rows];

	for (int i = 0; i < columns; i++) {
		for (int j = 0; j < columns; j++) {
			getColumn(dataMatrix, rows, i, column1);
			getColumn(dataMatrix, rows, j, column2);
			matrix[i][j] = computeCorrelation(column1, column2, rows);
		}
	}
	
	return matrix;
}

/************* Additional tests **************/

void runStudentTests6 (double **matrix, int rows, int columns) 
{
    std::cout << "--------------------------------------------------------------------------------" << std::endl
              << "Your own tests could print something here. See the runStudentTests6() function" << std::endl;

    // You can insert any code here — it will be called in the demo mode of Exercise 6

    std::cout << "--------------------------------------------------------------------------------" << std::endl;
}
