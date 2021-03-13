#pragma once
#include <string>
#include <vector>

bool testFunction (double (*func)(double[], std::vector<int>&),
                   std::string bad[], int nbBad, 
                   std::string good[], int nbGood,
                   std::string &path, double eps, bool print);
