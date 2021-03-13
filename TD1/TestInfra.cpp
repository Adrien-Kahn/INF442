
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <cmath>

/**************************** Tests - DO NOT TOUCH! ****************************/

std::ostream &printArray (std::ostream &str, double values[], int length, std::string prefix = "") {
  str << prefix << "Length: " << length << std::endl;

  str << prefix;
  for (int i = 0; i < length; i++)
    str << values[i] << " ";

  return str << std::endl;
}

void loadScenarii (std::ifstream &is, int nb_scenarii, int size_scenario, 
                  std::vector<std::vector<int> > &scenarii, std::vector<double> &results) {
  assert(is.is_open());

  for (int i = 0; i < nb_scenarii; i++) {
    for (int j = 0; j < size_scenario; j++) 
      is >> scenarii[i][j];

    is >> results[i];
  }
} 

void loadArray (std::ifstream &is, int length, int max, double values[]) {
  assert(is.is_open());
  for (int i = 0; i < length; i++) {
    is >> values[i]; 
  }
  for (int i = length; i < max; i++) {
    values[i] = 0;
  }
}

bool loadTests(std::string filename, int &nb_scenarii, int &size_scenario,
               std::vector<std::vector<int> > &scenarii,
               std::vector<double> &results, int &length, double *&values) {
  std::ifstream is(filename);
  if (!is.is_open())
    return false;

  int shift_index;
  is >> nb_scenarii >> size_scenario >> shift_index;

  scenarii.clear();
  scenarii.resize(nb_scenarii);
  for (int i = 0; i < nb_scenarii; i++)
    scenarii[i].resize(size_scenario);

  results.clear();
  results.resize(nb_scenarii);
  loadScenarii(is, nb_scenarii, size_scenario, scenarii, results);

  int maxShift = 0;
  if (shift_index >= 0) 
    for (int i = 0; i < nb_scenarii; i++)
      if (scenarii[i][shift_index] > maxShift)
        maxShift = scenarii[i][shift_index];

  is >> length;
  values = new double [length + maxShift];
  loadArray(is, length, length + maxShift, values);

  is.close();
  return true;
}

const int lengthIndex = 0;

bool testScenarii (double (*func)(double[], std::vector<int>&), std::string &filename, double eps, bool print) {
  int nb_scenarii;
  int size_scenario;
  std::vector<double> results;
  int length;
  double *values;
  std::vector<std::vector<int> > scenarii;
  
  if (!loadTests(filename, nb_scenarii, size_scenario, scenarii, results, length, values)) {
    std::cout << "\tProblem loading tests from the file " << filename << std::endl;
    return false;
  }

  if (print) {
    std::cout << "\tRead " << nb_scenarii << " scenarii of size " << size_scenario 
              << " from the file " << filename << std::endl
              << "\tWill test on the following array" << std::endl;
    printArray(std::cout, values, length, "\t");
  }

  bool found = false;
  for (int i = 0; i < nb_scenarii; i++) {
    double stat = func(values, scenarii[i]);
    // std::cout << stat << " " << fabs(stat - results[i]) << std::endl;
    if (fabs(stat - results[i]) < eps) {
      if (print) {
        std::cout << std::endl << "\tFound scenario ( ";
        for (int j = 0; j < size_scenario; j++)
          std::cout << scenarii[i][j] << " ";

        std::cout << ") matching the computed value " << stat
                  << std::endl;
      }
      
      found = true;
      break;
    }
  }

  delete[] values;
  return found;
}

bool testScenariiMsg (double (*func)(double[], std::vector<int>&), std::string scenarii[], int nb, bool needAll, 
                      std::string &path, double eps, bool print) {
  bool found = false;
  bool all = true;

  if (print)
    std::cout << std::endl << "\tWill test on " << (nb/2) << " scenario file(s) now..." << std::endl;

  // If one found, keep looking only if need to satisfy all
  for (int i = 0; (!found || needAll) && i < nb; i += 2) {
    std::string filename = path + scenarii[i];
    std::string message = scenarii[i+1];

    if (print)
      std::cout << std::endl;

    found = testScenarii(func, filename, eps, print);
    all = all && found;
    if (found && print && !message.empty())
        std::cerr << "\t" << message << std::endl;
  }

  return needAll ? all : found;
}

bool testFunction (double (*func)(double[], std::vector<int>&),
                   std::string bad[], int nbBad, 
                   std::string good[], int nbGood,
                   std::string &path, double eps, bool print) 
{
  if (print)
    std::cout << std::endl;

  // For the sake of clarity in the calls below
  bool needAll = true; 
  bool needOne = false;

  if (testScenariiMsg(func, bad, nbBad, needOne, path, eps, print)) {
    return false;
  }

  if (testScenariiMsg(func, good, nbGood, needAll, path, eps, print)) {
    return true;
  }

  std::cout << "\tYou have either made a mistake we haven't thought of or several mistakes." << std::endl;
  return false;
}
