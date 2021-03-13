#include <vector>

double hash (double vector[], int length) 
{
    double sum  = 0;
    for (int i = 0; i < length; i++)
        sum += vector[i];
    return sum;
} 

double **unpack (double values[], std::vector<int> &scenario, int size, int rows, int columns) 
{
    double **matrix = new double*[size];
    for (int i = 0; i < size; i++) {
        matrix[i] = new double[size];
        if (i < rows) {
            for (int j = 0; j < columns; j++)
                matrix[i][j] = values[i*columns + j];
            for (int j = columns; j < size; j++)
                matrix[i][j] = 1;
        } else {
            for (int j = 0; j < size; j++)
                matrix[i][j] = 1;
        }
    }

    return matrix;
}

void cleanup (double **matrix, int rows)
{
    for (int i = 0; i < rows; i++)
        delete[] matrix[i];
    delete matrix;
}

double **matrixFromArray (double array[], int rows, int columns) 
{
    double **matrix = new double*[rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[columns];
        for (int j = 0; j < columns; j++)
            matrix[i][j] = array[i*columns + j];
    }
    return matrix;
}
