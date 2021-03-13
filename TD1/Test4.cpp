#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include "TestInfra.hpp"
#include "Stats.hpp"
#include "ArrayToMatrix.hpp"

const bool defaultPrint = true;
const std::string defaultPath = "./tests/";
const double defaultEps = 0.001;

int findSize (double values[], int length, double &eps) {
    for (int i = length - 1; i >= 0; i--)
        if (fabs(values[i]) > eps)
            return i + 1;

    return 0;
}

double mockGetRow (double values[], std::vector<int> &scenario) {
    int rows = scenario[0];
    int columns = scenario[1];
    int size = scenario[2];
    int dStart = scenario[3];
    int dEnd = scenario[4];
    double eps = 1.0/pow(10, scenario[5]);

    double **matrix = unpack(values, scenario, size, rows, columns);
    double row[size];

    double mHash = 1;
    for (int i = 0; i < rows; i++) {
        memset(row, 0, size*sizeof(double));
        getRow(matrix, columns, i, row);
        mHash *= hash(row, findSize(row, size, eps));
    }

    cleanup(matrix, size);
        
    return mHash;
}

double mockGetColumn (double values[], std::vector<int> &scenario) {
    int rows = scenario[0];
    int columns = scenario[1];
    int size = scenario[2];
    int dStart = scenario[3];
    int dEnd = scenario[4];
    double eps = 1.0/pow(10, scenario[5]);

    double **matrix = unpack(values, scenario, size, rows, columns);
    double column[size];

    double mHash = 1;
    for (int i = 0; i < columns; i++) {
        memset(column, 0, size*sizeof(double));
        getColumn(matrix, rows, i, column);
        mHash *= hash(column, findSize(column, size, eps));
    }

    cleanup(matrix, size);
        
    return mHash;
}

bool testGetRow(std::string path, double eps, bool print)
{
  std::cout << "Testing the function getRow()...\t\t\t\t";

  std::string bad[] = {
    "getRowShift.txt", "Check the loop bounds: the indices of an array of size n go from 0 to n-1!",
    "getRowSwitch.txt", "Check the order of indices: did you invert columns and rows?",
    "getRowErase.txt", "Did you use the 'index' parameter on the wrong side of the assignment?",
};
  int nbBad = sizeof(bad)/sizeof(bad[0]);

  std::string good[] = {
    "getRow0.txt", ""
  };
  int nbGood = sizeof(good)/sizeof(good[0]);

  bool result = testFunction(mockGetRow, bad, nbBad, good, nbGood, path, eps, print);
  std::cout << (result ? "[OK]" : "[NOK]") << std::endl;
  return result;
}

bool testGetColumn(std::string path, double eps, bool print)
{
  std::cout << "Testing the function getColumn()...\t\t\t\t";

  std::string bad[] = {
    "getColumnShift.txt", "Check the loop bounds: the indices of an array of size n go from 0 to n-1!",
    "getColumnSwitch.txt", "Check the order of indices: did you invert columns and rows?",
    "getColumnErase.txt", "Did you use the 'index' parameter on the wrong side of the assignment?",
  };
  int nbBad = sizeof(bad)/sizeof(bad[0]);

  std::string good[] = {
    "getColumn0.txt", ""
  };
  int nbGood = sizeof(good)/sizeof(good[0]);

  bool result = testFunction(mockGetColumn, bad, nbBad, good, nbGood, path, eps, print);
  std::cout << (result ? "[OK]" : "[NOK]") << std::endl;
  return result;
}

void demoMatrix (double **matrix, int rows, int columns)
{
    std::cout << "Matrix:" << std::endl;
    printMatrix(matrix, rows, columns);
    std::cout << std::endl;

    double row[columns];
    double tempR[columns];
    memset(tempR, 0, sizeof(tempR));
    for (int i = 0; i < rows; i++) {
        memset(row, 0, sizeof(row));
        getRow(matrix, columns, i, row);
        if (memcmp(row, tempR, sizeof(row)) != 0)
        {
            std::cout << "Row " << i << ": ";
            printArray(row, columns);
            std::cout << std::endl;
        }
    }

    std::cout << std::endl;

    double column[rows];
    double tempC[rows];
    memset(tempC, 0, sizeof(tempC));
    for (int i = 0; i < columns; i++) {
        memset(column, 0, sizeof(column));
        getColumn(matrix, rows, i, column);
        if (memcmp(column, tempC, sizeof(column)) != 0)
        {
            std::cout << "Column " << i << ": ";
            printArray(column, rows);
            std::cout << std::endl;
        }
    }
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

    if (!tests) 
    {
        // If the program is executed without additional arguments do a demo
        std::cout << std::endl << "Running in demo mode" << std::endl;
        printUsage(argv[0]);
        std::cout << std::endl;

        demoMatrix(matrix, rows, columns);
        runStudentTests4(matrix, rows, columns);

        return 0;
    }
    else
    {
        // If additional arguments are provided run the tests
        std::cout << std::endl << "Running in the test mode" << std::endl;

        return !testGetRow(path, eps, print) ? 1 :
               !testGetColumn(path, eps, print) ? 2 : 
               0;
    }

    cleanup(matrix, rows);
}
