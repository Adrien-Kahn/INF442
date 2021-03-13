#include <iostream>
#include <iomanip> 
#include <cmath>
#include <string>

#include "Stats.hpp"

using namespace::std;

double **transpose (double **matrix, int rows, int columns) {
    double **transposed = new double*[columns];
    for (int i = 0; i < columns; i++) {
        transposed[i] = new double[rows];
        getColumn(matrix, rows, i, transposed[i]);
    }
    return transposed;
}

void cleanup (double **matrix, int rows)
{
    for (int i = 0; i < rows; i++)
        delete[] matrix[i];
    delete matrix;
}

int main () {
    int rows, columns;
    string r_name, c_name;

    // Read in the matrix dimensions
    cin >> columns >> c_name;
    cin >> rows >> r_name;
    
    // Reserve memory for the rows (each is a 'double*' pointer)
    double **matrix = new double*[rows];
    // Read in the matrix
    readMatrix(matrix, rows, columns);

    // Confirm
    cout << "Read in the "
         << rows << " " << r_name << " x " 
         << columns << " " << c_name 
         << " matrix" << endl
         ;
    printMatrix(matrix, rows, columns);

    cout << endl;

    // Print a table of means and variances of the rows
    cout << r_name << endl 
         << "\tmean\tvariance\tsample variance" << endl
         << "\t-------\t------------\t---------------" << endl
         ;
    for (int i = 0; i < rows; i++) {
        double row[columns];
        getRow(matrix, columns, i, row);
        cout << setw(2) << fixed << setprecision(1)
             << i << "\t"
             << computeMean(row, columns) << "\t"
             << computeVariance(row, columns) << "\t"
             << computeSampleVariance(row, columns) 
             << endl
             ;
    }

    cout << endl;

    // Print a table of means and variances of the columns
    cout << c_name << endl 
         << "\tmean\tvariance\tsample variance" << endl
         << "\t-------\t------------\t---------------" << endl
         ;
    for (int i = 0; i < columns; i++) {
        double column[rows];
        getColumn(matrix, rows, i, column);
        cout << setw(2) << fixed << setprecision(1)
             << i << "\t"
             << computeMean(column, rows) << "\t"
             << computeVariance(column, rows) << "\t"
             << computeSampleVariance(column, rows) 
             << endl
             ;
    }

    cout << endl;

    // Compute the column correlation matrix
    double **correlation = computeCorrelationMatrix(matrix, rows, columns);
    // Correlation matrices are square
    cout << c_name << "' correlation matrix (" << columns << " by " << columns << ")" << endl;
    printMatrix(correlation, columns, columns);
    // Free the memory occupied by the correlation matrix
    cleanup(correlation, columns);

    cout << endl;

    // Transpose the matrix
    double **transposed = transpose(matrix, rows, columns);
    // Compute the row correlation matrix (column correlation of the transposed matrix)
    correlation = computeCorrelationMatrix(transposed, columns, rows);
    // Correlation matrices are square
    cout << r_name << "' correlation matrix (" << rows << " by " << rows << ")" << endl;
    printMatrix(correlation, rows, rows);
    // Free the memory occupied by the transposed matrix
    cleanup(transposed, columns);
    // Free the memory occupied by the correlation matrix
    cleanup(correlation, rows);

    return 0;
}
