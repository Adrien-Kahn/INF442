#pragma once

#include "dataset.hpp"
#include "neuron.hpp"

const int default_nb_neurons = 5;
const int default_nb_epochs = 100;
const double default_learning_rate = 0.1;

class OneLayerPerceptron
{
protected:
    int dim;
    int size;
    std::vector<Node *> inputs;
    std::vector<Neuron *> hidden;
    Neuron *output;

    int epoch;
    double rate;
    double decay;

public:
    OneLayerPerceptron(int _dim, int _size, double _rate, double _decay,
                       std::function<double(double)> _activation,
                       std::function<double(double)> _activation_der);
    ~OneLayerPerceptron();

    int getNbNeurons();

    double getLearningRate();
    double getDecay();

    // "Using Learning Rate Schedules for Deep Learning Models
    // in Python with Keras", by Jason Brownlee
    // https://machinelearningmastery.com/using-learning-rate-schedules-deep-learning-models-python-keras/
    void setLearningRate(double _rate);
    void initLearningRate(double _rate);
    void decayLearningRate();

protected:
    double normalise(double val, Dataset *data, int coord);
    double denormalise(double val, Dataset *data, int coord);

    virtual void prepareInputs(Dataset *data, int row, int regr, bool print = false);
    virtual void computeHiddenStep(bool print = false);
    virtual double computeOutputStep(Dataset *data, int row, int regr, bool print = false);
    virtual void propagateBackHidden(bool print = false);

public:
    virtual double run(Dataset *data, int row, int regr, bool print = false);
};
