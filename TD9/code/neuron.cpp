#include "neuron.hpp"
#include "node.hpp"

#include <random>
#include <chrono> // To seed the random generator
#include <cmath>
#include <cassert>
#include <iostream>

// Vanilla constructor
Neuron::Neuron(int _nb_dendrites)
    : nb_dendrites(_nb_dendrites),
      dendrites(_nb_dendrites + 1, nullptr),
      axon(new Node()),
      activation([](double x) -> double { return (x > 0) ? 1 : 0; }),
      activation_der([](double x) -> double { return 0; }),
      weights(_nb_dendrites + 1, 1)
{
    // Setup "bias-related elements"
    dendrites[0] = new Node(-1);
    weights[0] = 0;
    initWeights();
    count++;
}

// Dendrites + activation function
Neuron::Neuron(int _nb_dendrites,
               std::vector<Node *> &_dendrites,
               std::function<double(double)> _activation,
               std::function<double(double)> _activation_der)
    : nb_dendrites(_nb_dendrites),
      dendrites(_dendrites),
      axon(new Node()),
      activation(_activation),
      activation_der(_activation_der),
      weights(_nb_dendrites + 1, 1)
{
    // Setup "bias-related elements"
    dendrites.insert(dendrites.begin(), new Node(-1));
    weights[0] = 0;
    initWeights();
    count++;
}

// Dendrites + activation function
Neuron::Neuron(int _nb_dendrites,
               std::vector<Neuron *> &_neurons,
               std::function<double(double)> _activation,
               std::function<double(double)> _activation_der)
    : nb_dendrites(_nb_dendrites),
      dendrites(_nb_dendrites + 1, nullptr),
      axon(new Node()),
      activation(_activation),
      activation_der(_activation_der),
      weights(_nb_dendrites + 1, 1)
{
    // Setup "bias-related elements"
    dendrites[0] = new Node(-1);
    weights[0] = 0;

    initWeights();

    for (int i = 0; i < nb_dendrites; i++)
        dendrites[i + 1] = _neurons[i]->axon;

    count++;
}

Neuron::~Neuron()
{
    delete dendrites[0];
    delete axon;
    count--;
}

void Neuron::initWeights()
{
    // We randomly initialise the weights for the dendrites,
    // not for the bias, the weight for the bias is initialised to 0
    // https://stats.stackexchange.com/questions/47590/what-are-good-initial-weights-in-a-neural-network
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();

    double upper_bound = 1.0 / std::sqrt(nb_dendrites);
    double lower_bound = -upper_bound;
    std::uniform_real_distribution<double>
        unif(lower_bound, upper_bound);
    std::default_random_engine re;

    // obtain a seed from the timer
    myclock::duration d = myclock::now() - beginning;
    re.seed(d.count());

    for (int i = 1; i <= nb_dendrites; i++)
        weights[i] = unif(re);
}

// Getters
int Neuron::getNbDendrites() const
{
    return nb_dendrites;
}

const Node *Neuron::getDendrite(int pos) const
{
    assert(0 <= pos && pos < nb_dendrites);
    return dendrites[pos + 1];
}

const Node *Neuron::getAxon() const
{
    return axon;
}

double Neuron::getWeight(int pos) const
{
    assert(0 <= pos && pos < nb_dendrites);
    return weights[pos + 1];
}

double Neuron::getBias() const
{
    return weights[0];
}

const Node *Neuron::getBiasDendrite() const
{
    return dendrites[0];
}

double Neuron::getLearningRate() const
{
    return rate;
}

double Neuron::getOutputSignal() const 
{
    return axon->getSignal();
}

double Neuron::getCollectedInput() const
{
    return collected_input;
}

// Setters
void Neuron::setWeights(std::vector<double> &_weights)
{
    assert(_weights.size() == nb_dendrites + 1);
    weights = _weights;
}

void Neuron::setWeight(int pos, double _weight)
{
    assert(0 <= pos && pos < nb_dendrites);
    weights[pos + 1] = _weight;
}

void Neuron::setBias(double _bias)
{
    weights[0] = _bias;
}

void Neuron::setLearningRate(double _rate)
{
    rate = _rate;
}

void Neuron::setCollectedInput(double _collected_input) {
    collected_input = _collected_input;
}

void Neuron::setBackValue(double back_value)
{
    axon->setBackValue(back_value);
}

// Forward propagation
void Neuron::step()
{

	// Checking that all nodes are connected
	for (int k = 0; k <= nb_dendrites; k++) {
		assert(dendrites[k]);	
	}
	assert(axon);

	double z = 0;
	
	for (int k = 0; k <= nb_dendrites; k++) {
		z += (dendrites[k])->getSignal() * weights[k];
	}
	
	double y = activation(z);
	
	collected_input = z;
	axon->setSignal(y);

}

// Backward propagation (training)
void Neuron::step_back()
{

	double err = axon->getBackValue() * activation_der(collected_input);
	
	for (int k = 1; k <= nb_dendrites; k++) {
		(dendrites[k])->setBackValue(err * weights[k]);
	}
	
	for (int k = 0; k <= nb_dendrites; k++) {
		weights[k] -= rate * err * dendrites[k]->getSignal();
	}

}

int Neuron::count = 0;

int Neuron::getCount() {
    return count;
}

/********************************************************************/
/* Auxiliary function for printing */

std::ostream &operator<<(std::ostream &str, const Neuron &n)
{
    using namespace std;

    int count = n.getNbDendrites();

    str << "Weights:";
    for (int i = 0; i < count; i++)
        str << " " << n.getWeight(i);
    str << endl;

    str << "Bias: " << n.getBias() << endl;

    str << "Dendrites:";
    for (int i = 0; i < count; i++)
    {
        if (n.getDendrite(i) != nullptr)
            str << " " << *(n.getDendrite(i));
        else
            str << " - ";
    }
    str << endl;

    str << "Bias: " << *n.getBiasDendrite() << endl;
    str << "Axon: " << *n.getAxon() << endl;

    return str;
}
