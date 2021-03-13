#include <iostream>
#include <string>
#include <vector>
#include "TestInfra.hpp"
#include "Stats.hpp"

const bool defaultPrint = true;
const std::string defaultPath = "./tests/";
const double defaultEps = 0.001;

double mockComputeMean (double values[], std::vector<int> &scenario) {
  return computeMean(values, scenario[0]);
}

bool testComputeMean (std::string path, double eps, bool print) {
  std::cout << "Testing the function computeMean()...\t\t\t\t";

  std::string bad[] = {
    "computeMeanShift.txt", "Check the loop bounds: the indices of an array of size n go from 0 to n-1!",
    "computeMeanLength.txt", "Check the divisor: is it equal to the number of values?",
    "computeMeanOW.txt", "Check the accumulator in the loop: did you add values or did you overwrite it?"
  };
  int nbBad = sizeof(bad)/sizeof(bad[0]);

  std::string good[] = {
    "computeMean0.txt", ""
  };
  int nbGood = sizeof(good)/sizeof(good[0]);

  bool result = testFunction(mockComputeMean, bad, nbBad, good, nbGood, path, eps, print);
  std::cout << (result ? "[OK]" : "[NOK]") << std::endl;
  return result;
}

double mockComputeVariance (double values[], std::vector<int> &scenario) {
  return computeVariance(values, scenario[0]);
}

bool testComputeVariance (std::string path, double eps, bool print) {
  std::cout << "Testing the function computeVariance()...\t\t\t";

  std::string bad[] = {
    "computeVarianceShift1.txt", "Check the loop bounds: the indices of an array of size n go from 0 to n-1!",
    "computeVarianceShift2.txt", "Check the loop bounds: the indices of an array of size n go from 0 to n-1!",
    "computeVarianceLength1.txt", "Check the divisor: is it equal to the number of values?",
    "computeVarianceLength2.txt", "Check the divisor: is it equal to the number of values?",
    "computeVarianceOW1.txt", "Check the accumulator in the loop: did you add values or did you overwrite it?",
    "computeVarianceOW2.txt", "Check the accumulator in the loop: did you add values or did you overwrite it?"
  };
  int nbBad = sizeof(bad)/sizeof(bad[0]);

  std::string good[] = {
    "computeVariance01.txt", "",
    "computeVariance02.txt", ""
  };
  int nbGood = sizeof(good)/sizeof(good[0]);

  bool result = testFunction(mockComputeVariance, bad, nbBad, good, nbGood, path, eps, print);
  std::cout << (result ? "[OK]" : "[NOK]") << std::endl;
  return result;
}

double mockComputeSampleVariance (double values[], std::vector<int> &scenario) {
  return computeSampleVariance(values, scenario[0]);
}

bool testComputeSampleVariance (std::string path, double eps, bool print) {
  std::cout << "Testing the function computeSampleVariance()...\t\t\t";

  std::string bad[] = {
    "computeSampleVariance1.txt", "Check the denominator in the definition of the (unbiased) sample variance",
  };
  int nbBad = sizeof(bad)/sizeof(bad[0]);

  std::string good[] = {
    "computeSampleVariance0.txt", ""
  };
  int nbGood = sizeof(good)/sizeof(good[0]);

  bool result = testFunction(mockComputeSampleVariance, bad, nbBad, good, nbGood, path, eps, print);
  std::cout << (result ? "[OK]" : "[NOK]") << std::endl;
  return result;
}

double mockComputeStandardDeviation (double values[], std::vector<int> &scenario) {
  return computeStandardDeviation(values, scenario[0]);
}

bool testComputeStandardDeviation (std::string path, double eps, bool print) {
  std::cout << "Testing the function computeStandardDeviation()...\t\t";

  std::string bad[] = {
    "computeStandardDeviation1.txt", "Did you take the square root correctly?",
    "computeStandardDeviation2.txt", "Did you use pow(..., 1/2) to take the square root? If so, recall that integer division gives 1/2 = 0.",
  };
  int nbBad = sizeof(bad)/sizeof(bad[0]);

  std::string good[] = {
    "computeStandardDeviation0.txt", ""
  };
  int nbGood = sizeof(good)/sizeof(good[0]);

  bool result = testFunction(mockComputeStandardDeviation, bad, nbBad, good, nbGood, path, eps, print);
  std::cout << (result ? "[OK]" : "[NOK]") << std::endl;
  return result;
}

double mockComputeSampleStandardDeviation (double values[], std::vector<int> &scenario) {
  return computeSampleStandardDeviation(values, scenario[0]);
}

bool testComputeSampleStandardDeviation (std::string path, double eps, bool print) {
  std::cout << "Testing the function computeSampleStandardDeviation()...\t";

  std::string bad[] = {
    "computeSampleStandardDeviation1.txt", "Did you take the square root correctly?",
    "computeStandardDeviation2.txt", "Did you use pow(..., 1/2) to take the square root? If so, recall that integer division gives 1/2 = 0.",
  };
  int nbBad = sizeof(bad)/sizeof(bad[0]);

  std::string good[] = {
    "computeSampleStandardDeviation0.txt", ""
  };
  int nbGood = sizeof(good)/sizeof(good[0]);

  bool result = testFunction(mockComputeSampleStandardDeviation, bad, nbBad, good, nbGood, path, eps, print);
  std::cout << (result ? "[OK]" : "[NOK]") << std::endl;
  return result;
}

void demoArrayStats (double sample[], int length)
{
    double mean = computeMean(sample, length);
    double variance = computeVariance(sample, length);
    double samplevariance = computeSampleVariance(sample, length);
    double sampledevi = computeSampleStandardDeviation(sample, length);
    double devi = computeStandardDeviation(sample, length);

    std::cout << "Sample: ";
    printArray(sample, length);
    std::cout << std::endl;

    std::cout << "Length\t\t\t\t: " << length << std::endl
              << "Mean\t\t\t\t: " << mean << std::endl
              << "Population variance\t\t: " << variance << std::endl
              << "Sample variance\t\t\t: " << samplevariance << std::endl
              << "Population standard deviation\t: " << devi << std::endl
              << "Sample standard deviation\t: " << sampledevi << std::endl
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

    if (!tests) 
    {
        // If the program is executed without additional arguments do a demo
        std::cout << std::endl << "Running in demo mode" << std::endl;
        printUsage(argv[0]);
        std::cout << std::endl;

        demoArrayStats(sample, sample_length);
        runStudentTests3(sample, sample_length);

        return 0;
    }
    else
    {
        // If additional arguments are provided run the tests
        if (print)
            std::cout << std::endl << "Running in the test mode" << std::endl;

        return !testComputeMean(path, eps, print) ? 1 : 
               !testComputeVariance(path, eps, print) ? 2 :
               !testComputeSampleVariance(path, eps, print) ? 3 : 
               !testComputeStandardDeviation(path, eps, print) ? 4 :
               !testComputeSampleStandardDeviation(path, eps, print) ? 5 :
               0;
    }
}
