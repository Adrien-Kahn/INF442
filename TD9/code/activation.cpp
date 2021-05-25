#include "activation.hpp"

double sigma(double x)
{
    return 1.0 / (1.0 + exp(-x));
}

double sigma_der(double x)
{
    double s = sigma(x);
    return s * (1 - s);
}

double id(double x)
{
    return x;
}

double id_der(double x)
{
    return 1;
}
