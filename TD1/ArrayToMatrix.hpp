#pragma once

#include <vector>

double hash (double vector[], int length);

double **unpack (double values[], std::vector<int> &scenario, int size, int rows, int columns);

void cleanup (double **matrix, int rows);

double **matrixFromArray (double array[], int rows, int columns);
