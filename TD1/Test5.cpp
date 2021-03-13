#include <iostream>
#include <string>
#include "TestInfra.hpp"
#include "Stats.hpp"
#include "ArrayToMatrix.hpp"

const bool defaultPrint = true;
const std::string defaultPath = "./tests/";
const double defaultEps = 0.001;

double mockComputeCovariance (double values[], std::vector<int> &scenario) {
  int length = scenario[0];
  return computeCovariance(values, values + length, length);
}

double mockComputeCorrelation (double values[], std::vector<int> &scenario) {
  int length = scenario[0];
  return computeCorrelation(values, values + length, length);
}

bool testComputeCovariance(std::string path, double eps, bool print)
{
  std::cout << "Testing the function computeCovariance()...\t";

  std::string bad[] = {
    "computeCovarianceShift.txt", "Check the loop bounds: the indices of an array of size n go from 0 to n-1!",
    "computeCovarianceLength.txt", "Check the divisor: is it equal to the number of values?",
    "computeCovarianceOW.txt", "Check the accumulator in the loop: did you add values or did you overwrite it?",
};
  int nbBad = sizeof(bad)/sizeof(bad[0]);

  std::string good[] = {
    "computeCovariance0.txt", ""
  };
  int nbGood = sizeof(good)/sizeof(good[0]);

  bool result = testFunction(mockComputeCovariance, bad, nbBad, good, nbGood, path, eps, print);
  std::cout << (result ? "[OK]" : "[NOK]") << std::endl;
  return result;
}

bool testComputeCorrelation(std::string path, double eps, bool print)
{
  std::cout << "Testing the function computeCorrelation()...\t";

  std::string bad[] = {
    "computeCorrelation1.txt", "Did you devide by both standard deviations? Check the / and * operations and parentheses.",
  };
  int nbBad = sizeof(bad)/sizeof(bad[0]);

  std::string good[] = {
    "computeCorrelation0.txt", ""
  };
  int nbGood = sizeof(good)/sizeof(good[0]);

  bool result = testFunction(mockComputeCorrelation, bad, nbBad, good, nbGood, path, eps, print);
  std::cout << (result ? "[OK]" : "[NOK]") << std::endl;
  return result;
}

void demoCovarianceCorrelation (double sample1[], double sample2[], int length)
{
    std::cout << "Sample 1: ";
    printArray(sample1, length);
    std::cout << std::endl
              << "Length:\t\t" << length << std::endl
              << "Mean:\t\t" << computeMean(sample1, length) << std::endl
              << "Variance:\t" << computeVariance(sample1, length) << std::endl
              << "Std dev:\t" << computeStandardDeviation(sample1, length) << std::endl
              << std::endl
              ;
              
    std::cout << "Sample 2: ";
    printArray(sample2, length);
    std::cout << std::endl
              << "Length:\t\t" << length << std::endl
              << "Mean:\t\t" << computeMean(sample2, length) << std::endl
              << "Variance:\t" << computeVariance(sample2, length) << std::endl
              << "Std dev:\t" << computeStandardDeviation(sample2, length) << std::endl
              << std::endl
              ;

    std::cout << "Covariance 1 & 1: " << computeCovariance(sample1, sample1, length) << std::endl
              << "Covariance 1 & 2: " << computeCovariance(sample1, sample2, length) << std::endl
              << "Covariance 2 & 1: " << computeCovariance(sample2, sample1, length) << std::endl
              << "Covariance 2 & 2: " << computeCovariance(sample2, sample2, length) << std::endl
              << std::endl;

    std::cout << "Correlation 1 & 1: " << computeCorrelation(sample1, sample1, length) << std::endl
              << "Correlation 1 & 2: " << computeCorrelation(sample1, sample2, length) << std::endl
              << "Correlation 2 & 1: " << computeCorrelation(sample2, sample1, length) << std::endl
              << "Correlation 2 & 2: " << computeCorrelation(sample2, sample2, length) << std::endl
              ;
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
    double sample1[] = {
        15.17, 16.94, 14.94, 14.99, 13.77, 13.75,
        12.67, 12.14, 12.59, 12.48, 14.81, 14.29,
        12.74, 12.52, 11.65,
    };
    double sample2[] = {
        12.24, 11.42, 12.25,
        12.72, 11.64, 11.09, 11.22, 11.50, 11.36,
        11.84, 12.18, 11.04, 10.90, 11.80, 11.84,
    };
    // We have to pass the length of the array as an argument to the function,
    // since it is not "stored" in the array itself.
    //
    // The size of the array (in memory) is the size of one element
    // times the number of elements
    int sample_length = sizeof(sample1) / sizeof(sample1[0]);

    if (!tests) 
    {
        // If the program is executed without additional arguments do a demo
        std::cout << std::endl << "Running in demo mode" << std::endl;
        printUsage(argv[0]);
        std::cout << std::endl;

        demoCovarianceCorrelation(sample1, sample2, sample_length);
        runStudentTests5(sample1, sample2, sample_length);

        return 0;
    }
    else
    {
        // If additional arguments are provided run the tests
        std::cout << std::endl << "Running in the test mode" << std::endl;

        return !testComputeCovariance(path, eps, print) ? 1 :
               !testComputeCorrelation(path, eps, print) ? 2 : 
               0;
        }
}
