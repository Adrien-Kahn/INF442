#include <iostream>
#include <string>
#include <cmath>
#include "TestInfra.hpp"
#include "ArrayToMatrix.hpp"
#include "Stats.hpp"

const bool defaultPrint = true;
const std::string defaultPath = "./tests/";
const double defaultEps = 0.001;

int findSizeArray (double values[], int length, double &eps) {
    for (int i = length - 1; i >= 0; i--)
        if (fabs(values[i]) > eps)
            return i + 1;

    return 0;
}

int findSize (double **matrix, int side, double &eps) 
{
    int sizeH = 0;
    int sizeV = 0;
    
    for (int i = 0; i < side; i++)
        for (int j = side - 1; j >= 0; j--)
            if (fabs(matrix[i][j]) > eps) {
                sizeH = (sizeH > j + 1) ? sizeH : j + 1;
                break;
            }

    for (int i = 0; i < side; i++)
        for (int j = side - 1; j >= 0; j--)
            if (fabs(matrix[j][i]) > eps) {
                sizeV = (sizeV > j + 1) ? sizeV : j + 1;
                break;
            }

    return sizeH * (side + 1) + sizeV;
}

double mockComputeCovarianceMatrix (double values[], std::vector<int> &scenario) {
    int rows = scenario[0];
    int columns = scenario[1];
    int bound1 = scenario[2];
    int bound2 = scenario[3];
    double eps = 1.0/pow(10, scenario[7]);

    int maxBound = (bound1 > bound2) ? bound1 : bound2;
    int max = (rows > columns ? rows : columns);
    max = (max > maxBound) ? max : maxBound;

    double **matrix = unpack(values, scenario, max, rows, columns);

    for (int i = 0; i < max; i++)
      for (int j = 0; j < max; j++)
        if (i >= rows || j >= columns)
          matrix[i][j] = i + j;

    double **computed = computeCovarianceMatrix(matrix, rows, columns);

    // int temp = findSize(computed, columns, eps);
    // int trueRows = temp % (columns + 1);
    // int trueColumns = temp / (columns + 1);
    int trueRows = columns;
    int trueColumns = columns;

    double prod = 1;
    for (int i = 0; i < trueRows; i++) {
        double sum = 0;
        for (int j = 0; j < trueColumns; j++)
            sum += computed[i][j];
        if (fabs(sum) > eps)
            prod *= sum;
    }

    cleanup(matrix, max);
    cleanup(computed, columns);

    return prod;
}

double mockComputeCorrelationMatrix (double values[], std::vector<int> &scenario) {
    int rows = scenario[0];
    int columns = scenario[1];
    int bound1 = scenario[2];
    int bound2 = scenario[3];
    double eps = 1.0/pow(10, scenario[7]);

    int maxBound = (bound1 > bound2) ? bound1 : bound2;
    int max = (rows > columns ? rows : columns);
    max = (max > maxBound) ? max : maxBound;

    double **matrix = unpack(values, scenario, max, rows, columns);

    for (int i = 0; i < max; i++)
      for (int j = 0; j < max; j++)
        if (i >= rows || j >= columns)
          matrix[i][j] = i + j;

    double **computed = computeCorrelationMatrix(matrix, rows, columns);

    int temp = findSize(computed, columns, eps);
    int trueRows = temp % (columns + 1);
    int trueColumns = temp / (columns + 1);

    double prod = 1;
    for (int i = 0; i < trueRows; i++) {
        double sum = 0;
        for (int j = 0; j < trueColumns; j++)
            sum += computed[i][j];
        if (sum != 0)
            prod *= sum;
    }

    cleanup(matrix, max);
    cleanup(computed, columns);

    return prod;
}

bool testComputeCovarianceMatrix(std::string path, double eps, bool print)
{
  std::cout << "Testing the function computeCovarianceMatrix()...\t\t\t\t";

  std::string bad[] = {
    "computeCovarianceMatrixShift.txt", "Check the loop bounds: the indices of an array of size n go from 0 to n-1!",
    "computeCovarianceMatrixColumn.txt", "Do you pass the right length to getColumn()?",
    "computeCovarianceMatrixFlip.txt", "Did you forget to fill out the symmetric part of the covariance matrix?",
    "computeCovarianceMatrixDiagonal.txt", "Did you set the values at the diagonal to 1 instead of the variance?",
};
  int nbBad = sizeof(bad)/sizeof(bad[0]);

  std::string good[] = {
    "computeCovarianceMatrix0.txt", ""
  };
  int nbGood = sizeof(good)/sizeof(good[0]);

  bool result = testFunction(mockComputeCovarianceMatrix, bad, nbBad, good, nbGood, path, eps, print);
  std::cout << (result ? "[OK]" : "[NOK]") << std::endl;
  return result;
}

bool testComputeCorrelationMatrix(std::string path, double eps, bool print)
{
  std::cout << "Testing the function computeCorrelationMatrix()...\t\t\t\t";

  std::string bad[] = {
    "computeCorrelationMatrixShift.txt", "Check the loop bounds: the indices of an array of size n go from 0 to n-1!",
    "computeCorrelationMatrixColumn.txt", "Do you pass the right length to getColumn()?",
    "computeCorrelationMatrixFlip.txt", "Did you forget to fill out the symmetric part of the correlation matrix?",
    "computeCorrelationMatrixDiagonal.txt", "Did you set the values at the diagonal to variance instead of 1?",
  };
  int nbBad = sizeof(bad)/sizeof(bad[0]);

  std::string good[] = {
    "computeCorrelationMatrix0.txt", ""
  };
  int nbGood = sizeof(good)/sizeof(good[0]);

  bool result = testFunction(mockComputeCorrelationMatrix, bad, nbBad, good, nbGood, path, eps, print);
  std::cout << (result ? "[OK]" : "[NOK]") << std::endl;
  return result;
}

void demoMatrix (double **matrix, int rows, int columns)
{
    std::cout << "Matrix:" << std::endl;
    printMatrix(matrix, rows, columns);
    std::cout << std::endl;

    std::cout << "Covariance matrix:" << std::endl;
    double **covariance = computeCovarianceMatrix(matrix, rows, columns);
    printMatrix(covariance, columns, columns);
    std::cout << std::endl;

    std::cout << "Correlation matrix:" << std::endl;
    double **correlation = computeCorrelationMatrix(matrix, rows, columns);
    printMatrix(correlation, columns, columns);
    std::cout << std::endl;

    cleanup(covariance, columns);
    cleanup(correlation, columns);
}

void printUsage(char *name) {
  std::cout << std::endl
            << "USAGE: To test your code, run the program as follows:" << std::endl
            << std::endl
            << name << " test <print> <path> <epsilon>" << std::endl
            << std::endl
            << "The last three parameters are optional but have to be given in that order." << std::endl
            << "<print>\t\t - an integer - test information is printed only if it is positive," << std::endl
            << "\t\t   default value is " << defaultPrint << std::endl
            << std::endl
            << "<path> \t\t - the path to the test files, default path is '" << defaultPath << "'" << std::endl
            << std::endl
            << "<epsilon>\t - a double number, which defines the comparison sensitivity" << std::endl
            << "\t\t   for doubles, default is " << defaultEps << std::endl
            ;
}

/**
 * The main function
 * 
 * @param argc the number of arguments passed on the command line
 * @param argv the array of arguments — its length is argc
 * @return if executed successfully return 0, otherwise return an error code 
 */
int main(int argc, char *argv[])
{
    // Read in the arguments
    int arg = 1; 
    const bool tests = (argc > arg) ; arg++;
    const bool print = (argc > arg) ? std::stoi(argv[arg]) > 0 : defaultPrint; arg++;
    std::string path = (argc > arg) ? argv[arg] : defaultPath; arg++;
    const double eps = (argc > arg) ? std::stod(argv[arg]) : defaultEps;

    if (path[path.length() - 1] != '/')
        path += "/";

    // Sample array of data
    double sample[] = {
        15.17, 16.94, 14.94, 14.99, 13.77, 13.75,
        12.67, 12.14, 12.59, 12.48, 14.81, 14.29,
        12.74, 12.52, 11.65, 12.24, 11.42, 12.25,
        12.72, 11.64, 11.09, 11.22, 11.50, 11.36,
        11.84, 12.18, 11.04, 10.90, 11.80, 11.84,
    };
    // We have to pass the length of the array as an argument to the function,
    // since it is not "stored" in the array itself.
    //
    // The size of the array (in memory) is the size of one element
    // times the number of elements
    int sample_length = sizeof(sample) / sizeof(sample[0]);

    // The same sample but arranged as a 5x6 matrix
    int rows = 5;
    int columns = 6;
    double **matrix = matrixFromArray(sample, rows, columns);

    int result = 0;  
    if (!tests) 
    {
        // If the program is executed without additional arguments do a demo
        std::cout << std::endl << "Running in demo mode" << std::endl;
        printUsage(argv[0]);
        std::cout << std::endl;

        demoMatrix(matrix, rows, columns);
        runStudentTests6(matrix, rows, columns);
    }
    else
    {
        // If additional arguments are provided run the tests
        std::cout << std::endl << "Running in the test mode" << std::endl;

        result = !testComputeCovarianceMatrix(path, eps, print) ? 1 :
                 !testComputeCorrelationMatrix(path, eps, print) ? 2 : 
                 0;
    }

    cleanup(matrix, rows);
    return result;
}
