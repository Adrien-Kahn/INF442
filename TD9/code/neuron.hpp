#pragma once

#include "node.hpp"

#include <iostream>
#include <functional>
#include <vector>

class Neuron
{
private:
    int nb_dendrites;
    std::vector<Node *> dendrites;
    Node *axon;

    // by default use the step function
    std::function<double(double)> activation;
    std::function<double(double)> activation_der;

    // by default use identity weights and no bias
    // vector should have length size + 1: [0] = bias
    std::vector<double> weights;

    // Collected input
    // used for both forward and backward propagation
    double collected_input = 0;
    double rate = 1;

    // For testing only
    static int count;
public:
    // Vanilla constructor
    Neuron(int _nb_dendrites);

    // Dendrites + activation function
    Neuron(int _nb_dendrites,
           std::vector<Node *> &_dendrites,
           std::function<double(double)> _activation,
           std::function<double(double)> _activation_der);
    Neuron(int _nb_dendrites,
           std::vector<Neuron *> &_neurons,
           std::function<double(double)> _activation,
           std::function<double(double)> _activation_der);

    ~Neuron();

    void initWeights();

    // Getters
    int getNbDendrites() const;
    const Node *getDendrite(int pos) const;
    const Node *getAxon() const;
    const Node *getBiasDendrite() const;
    double getWeight(int pos) const;
    double getBias() const;
    double getLearningRate() const;
    double getOutputSignal() const;
    double getCollectedInput() const; // for testing only
    static int getCount(); // for testing only

    // Setters
    void setWeights(std::vector<double> &_weights);
    void setWeight(int pos, double _weight);
    void setBias(double _bias);
    void setLearningRate(double _rate);
    void setCollectedInput(double _collected_input); // for testing only

    void setBackValue(double back_value);

    // Forward propagation
    void step();
    // Backward propagation (training)
    void step_back();
};

std::ostream &operator<<(std::ostream &str, const Neuron &n);
