#pragma once

#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>

class Dataset
{
public:
    Dataset(std::ifstream &fin);
    ~Dataset();

    void Show(bool verbose);
    std::vector<double> &getInstance(int i);
    int getNbrSamples() const;
    int getDim() const;

    std::vector<double> *getMins();
    std::vector<double> *getMaxs();
    void mergeStats (Dataset &that);

    double getMin(int i);
    double getMax(int i);
    double getMean(int i);
    double getVariance(int i);
    double getStddev(int i);

private:
    int dim;
    int nsamples;
    std::vector<std::vector<double> > instances;

    std::vector<double> min;
    std::vector<double> max;
    std::vector<double> mean;
    std::vector<double> variance;
    std::vector<double> stddev;

    void computeStats();
};

std::ostream &operator<<(std::ostream &str, const std::vector<double> &vx);
