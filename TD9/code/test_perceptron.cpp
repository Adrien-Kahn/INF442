#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>

#include <assert.h>
#include "dataset.hpp"
#include "perceptron.hpp"
#include "activation.hpp"

using namespace std;

class OLP_test : public OneLayerPerceptron
{
public:
    enum OLPERR
    {
        OK,
        NO,
        SIZE,
        INPUTS,
        ACT,
        ACTDER
    };
    static double eps;

    OLP_test(int _dim, int _size, double _rate, double _decay) : OneLayerPerceptron(_dim, _size, _rate, _decay, sigma, sigma_der) {}
    OLP_test(int _dim, int _size, double _rate, double _decay, 
        std::function<double(double)> _activation, std::function<double(double)> _activation_der) : 
        OneLayerPerceptron(_dim, _size, _rate, _decay, _activation, _activation_der) {}

    bool test_inputs();
    OLPERR test_hidden();
    OLPERR test_output();
    bool test_act(std::function<double(double)> a, std::function<double(double)> da);

    static bool interpret_neuron_test(OLPERR code, std::string type);

    bool test_prepareInputs(Dataset *data, int row, int regr);
    bool test_computeHiddenStep();
    bool test_computeOutputStep(Dataset *data, int row, int regr);
    bool test_propagateBackHidden();
};

double OLP_test::eps = 0.0001;

bool OLP_test::test_inputs()
{
    for (int i = 0; i < dim; i++)
        if (inputs[i] == nullptr)
            return false;

    return true;
}

OLP_test::OLPERR OLP_test::test_hidden()
{
    for (int i = 0; i < size; i++)
        if (hidden[i] == nullptr)
            return NO;

    for (int i = 0; i < size; i++)
        if (hidden[i]->getNbDendrites() != dim)
            return SIZE;

    for (int i = 0; i < size; i++)
        for (int j = 0; j < dim; j++)
            if (hidden[i]->getDendrite(j) != inputs[j])
                return INPUTS;

    return OK;
}

OLP_test::OLPERR OLP_test::test_output()
{
    if (output == nullptr)
        return NO;

    if (output->getNbDendrites() != size)
        return SIZE;

    for (int i = 0; i < size; i++)
        if (output->getDendrite(i) != hidden[i]->getAxon())
            return INPUTS;

    return OK;
}

bool OLP_test::interpret_neuron_test(OLP_test::OLPERR code, std::string type)
{
    switch (code)
    {
    case OLP_test::NO:
        cout << "[NOK] - " << type << " neuron is not initialised correctly" << endl;
        return false;

    case OLP_test::SIZE:
        cout << "[NOK] - " << type << " neuron has incorrect number of dendrites" << endl;
        return false;

    case OLP_test::INPUTS:
        cout << "[NOK] - " << type << " neuron is not correctly connected to inputs" << endl;
        return false;

    case OLP_test::OK:
        return true;

    default:
        cout << "[NOK] - " << type << " neuron has an unknown problem" << endl;
        return false;
    }
}

bool OLP_test::test_act(std::function<double(double)> a, std::function<double(double)> da)
{
    for (int i = 0; i < dim; i++)
    {
        inputs[i]->setSignal(1);
        for (int j = 0; j < size; j++)
            hidden[j]->setWeight(i, 1);
    }
    double ad = a(dim);
    for (int i = 0; i < size; i++)
    {
        output->setWeight(i, 1);
        hidden[i]->step();
        if (fabs(hidden[i]->getAxon()->getSignal() - ad) > eps)
        {
            cout << "[NOK] - Incorrect activation function used for a hidden neuron (did you hardcode sigma instead of using the argument?)" << endl;
            return false;
        }
    }
    output->step();
    if (fabs(output->getAxon()->getSignal() - ad * size) > eps)
    {
        cout << "[NOK] - Incorrect activation function used for the output neuron (did you hardcode sigma instead of using the argument?)" << endl;
        return false;
    }
    output->setBackValue(1);
    output->step_back();
    for (int i = 0; i < size; i++)
    {
        if (fabs(output->getWeight(i) + ad - 1) > eps)
        {
            cout << "[NOK] - Incorrect activation function derivative used for the output neuron (did you hardcode sigma_der instead of using the argument?)" << endl;
            return false;
        }
        hidden[i]->setBackValue(1);
        hidden[i]->step_back();
    }
    double dad = da(dim);
    for (int i = 0; i < size; i++)
        for (int j = 0; j < dim; j++)
            if (fabs(hidden[i]->getWeight(j) + dad - 1) > eps)
            {
                cout << "[NOK] - Incorrect activation function derivative used for a hidden neuron (did you hardcode sigma_der instead of using the argument?)" << endl;
                return false;
            }
    return true;
}

bool OLP_test::test_prepareInputs(Dataset *data, int row, int regr) {
    prepareInputs(data, row, regr);

    std::vector<double> *mins = data->getMins();
    std::vector<double> *maxs = data->getMaxs();

    for (int j = 0; j < size; j++) 
        for (int i = 0; i < dim; i++) 
            hidden[j]->setWeight(i, 1);

    if (regr > 0) {
        double value = data->getInstance(row)[regr-1];
        double norm = normalise(value, data, regr-1);

        if (fabs(inputs[regr-1]->getSignal() - value) < eps) {
            cout << "[NOK] - Input data not normalised" << endl;
            return false;
        }

        if (fabs(inputs[regr-1]->getSignal() - norm) > eps) {
            cout << "[NOK] - Input data not initialised correctly" << endl;
            return false;
        }
    }

    if (regr <= dim) {
        double value = data->getInstance(row)[regr+1];
        double shifted_value = data->getInstance(row)[regr];
        double norm = normalise(value, data, regr + 1);
        double shifted_norm = normalise(shifted_value, data, regr);

        if (fabs(inputs[regr]->getSignal() - value) < eps) {
            cout << "[NOK] - Input data not normalised" << endl;
            return false;
        }

        if (fabs(inputs[regr]->getSignal() - shifted_value) < eps ||
                fabs(inputs[regr]->getSignal() - shifted_norm) < eps) {
            cout << "[NOK] - Forgot to skip the regression column" << endl;
            return false;
        } 

        if (fabs(inputs[regr]->getSignal() - norm) > eps) {
            cout << "[NOK] - Input data not initialised correctly" << endl;
            return false;
        }
    }

    return true;
}

bool OLP_test::test_computeHiddenStep() {
    double ax[size];
    for (int i = 0; i < size; i++)
        ax[i] = hidden[i]->getAxon()->getSignal();
    computeHiddenStep();
    bool result = true;
    for (int i = 0; i < size; i++) 
        result = result && hidden[i]->getAxon()->getSignal() != ax[i];
    return result;
}

bool OLP_test::test_computeOutputStep(Dataset *data, int row, int regr) {
    std::vector<double> *mins = data->getMins();
    std::vector<double> *maxs = data->getMaxs();
    double out = output->getAxon()->getSignal();
    double back = output->getAxon()->getBackValue();
    double den[size];
    for (int i = 0; i < size; i++)
        den[i] = output->getDendrite(i)->getBackValue();
    double output_sig = computeOutputStep(data, row, regr);
    if (fabs(output_sig - denormalise(output->getAxon()->getSignal(), data, regr) > eps)) {
        cout << "[NOK] - Incorrect output signal (probably forgot to set 'ret')" << endl;
        return false;
    }
    if (fabs(output->getAxon()->getSignal() - out) < eps) {
        cout << "[NOK] - Output signal not updated" << endl;
        return false;
    }
    if (fabs(output->getAxon()->getBackValue() - back) < eps) {
        cout << "[NOK] - Back propagation not initialised in the output neuron" << endl;
        return false;
    }
    if (fabs(output->getAxon()->getBackValue() - output->getAxon()->getSignal()) < eps) {
        cout << "[NOK] - Incorrectly initialised back propagation in the output neuron ('ret' instead of 'error'?)" << endl;
        return false;
    }
    for (int i = 0; i < size; i++) 
        if (hidden[i]->getAxon()->getBackValue() == den[i]) {
            cout << "[NOK] - Back propagation through the output neuron not done" << endl;
            return false;
        }

    return true;
}

bool OLP_test::test_propagateBackHidden() {
    double ax[dim];
    for (int i = 0; i < dim; i++)
        ax[i] = inputs[i]->getBackValue();
    for (int i = 0; i < size; i++)
        hidden[i]->setBackValue(100);
    propagateBackHidden();
    bool result = true;
    for (int i = 0; i < dim; i++) 
        if (fabs(inputs[i]->getBackValue() - ax[i]) < eps) {
            result = false;
            break;
        }
    return result;
}

void run_on_data(OneLayerPerceptron *perceptron,
                 Dataset &data,
                 int rounds,
                 int regr,
                 bool train = true,
                 bool print = false)
{
    clock_t tic, toc;
    clock_t cumulative = 0;
    double rss = 0;

    int rows = data.getNbrSamples();

    for (int round = 0; round < rounds; round++)
    {
        tic = clock();
        for (int row = 0; row < rows; row++)
        {
            if (print)
                cout << "Epoch/round = " << round
                     << ", row = " << row << endl;
            // double output = 0;
            double output = perceptron->run(&data, row, regr, print);
            if (print)
                cout << "\tOutput: " << output;
            double err = (output - data.getInstance(row)[regr]);
            if (print)
                cout << "\tError: " << err << endl;
            rss += err * err;
            if (print)
                cout << "Current RSS: " << rss << endl;
        }

        if (train)
        {
            if (print)
                cout << "Updating learning rate...\t";
            perceptron->decayLearningRate();
            if (print)
                cout << perceptron->getLearningRate() << endl;
        }

        toc = clock();
        cumulative += toc - tic;
    }

    cout << "Mean RSS: " << rss / rounds << endl;
    cout << "Total time elapsed = "
         << cumulative / ((float)(CLOCKS_PER_SEC))
         << "s" << endl;
    cout << "Mean time per epoch/round = "
         << cumulative / ((float)(CLOCKS_PER_SEC)) / rounds
         << "s" << endl;
}

bool test3()
{
    const int dim = 3;
    const int size = 5;
    const double rate = 1;
    const double decay = 0;

    cout << "Testing OneLayerPerceptron constructor...\t";
    OLP_test *perceptron = new OLP_test(dim, size, rate, decay, [](double x) -> double {return x*x;}, [](double x) -> double {return 2*x;});

    if (!perceptron->test_inputs())
    {
        cout << "[NOK] - Input nodes are not initialised correctly" << endl;
        return false;
    }

    if (!OLP_test::interpret_neuron_test(perceptron->test_hidden(), "A hidden"))
        return false;

    if (!OLP_test::interpret_neuron_test(perceptron->test_output(), "The output"))
        return false;

    if (!perceptron->test_act([](double x) -> double {return x*x;}, [](double x) -> double {return 2*x;}))
        return false;

    cout << "[OK]" << endl;

    cout << "Testing OneLayerPerceptron destructor...\t";

    delete perceptron;

    switch (Neuron::getCount())
    {
    case 1:
        cout << "[NOK] - Output neuron was not deleted" << endl;
        return false;

    case 0:
        break;

    default:
        cout << "[NOK] - Hidden neurons were not deleted" << endl;
        return false;
    }

    switch (Node::getCount())
    {
    case dim:
        cout << "[NOK] - Inputs were not deleted" << endl;
        return false;

    case 0:
        cout << "[OK]" << endl;
        break;

    default:
        cout << "[NOK] - Some nodes were not deleted" << endl;
        return false;
    }

    return true;
}

bool test4(const char *filename) {
    cout << "Preparing to test `OneLayerPerceptron::run`..." << endl;
    std::ifstream ftest(filename);
    if (ftest.fail())
    {
        std::cout << "\tCould not read from the file `" << filename << "` - make sure its available" 
                  << std::endl;
        return false;
    }

    Dataset test(ftest);
    const int dim = test.getDim() - 1;
    const int count = test.getNbrSamples();
    const int row = 1;
    cout << "\tRead test data from " << filename << endl;
    cout << "\t" << count << " rows of dimension " << test.getDim()
         << endl;

    const int regr = 2;
    const int size = default_nb_neurons;
    const double rate = 1;
    const double decay = 0;

    cout << "\tSetting up a perceptron for testing...\t";
    OLP_test perceptron(dim, size, rate, decay);
    cout << "done." << endl;

    cout << "...preparation done." << endl;

    cout << "Testing OneLayerPerceptron::prepareInputs()...\t\t";
    if (perceptron.test_prepareInputs(&test, row, regr)) 
        cout << "[OK]" << endl;
    else 
        return false;

    cout << "Testing OneLayerPerceptron::computeHiddenStep()...\t";
    if (!perceptron.test_computeHiddenStep()) {
        cout << "[NOK]" << endl;
        return false;
    }
    else
        cout << "[OK]" << endl;

    cout << "Testing OneLayerPerceptron::computeOutputStep()...\t";
    if (perceptron.test_computeOutputStep(&test, row, regr)) {
        cout << "[OK]" << endl;
    }
    else
        return false;

    cout << "Testing OneLayerPerceptron::propagateBackHidden()...\t";
    if (!perceptron.test_propagateBackHidden()) {
        cout << "[NOK]" << endl;
        return false;
    }
    else
        cout << "[OK]" << endl;

    return true;
}

int run(int argc, char *argv[])
{
    std::ifstream ftrain(argv[1]);
    std::ifstream ftest(argv[2]);

    if (ftrain.fail())
    {
        std::cout << "Cannot read from the testing file" << std::endl;
        return 1;
    }

    if (ftest.fail())
    {
        std::cout << "Cannot read from the testing file" << std::endl;
        return 1;
    }

    Dataset training(ftrain);
    int dim = training.getDim() - 1;
    int count = training.getNbrSamples();
    cout << "Read training data from " << argv[1] << endl;
    cout << count << " rows of dimension " << training.getDim()
         << endl;

    Dataset testing(ftest);
    cout << "Read testing data from " << argv[2] << endl;
    cout << count << " rows of dimension " << testing.getDim()
         << endl;

    assert(training.getDim() == testing.getDim());

    int regr = (argc > 3) ? std::atoi(argv[3]) : dim;
    int size = (argc > 4) ? std::atoi(argv[4]) : default_nb_neurons;
    int epochs = (argc > 5) ? std::atoi(argv[5]) : default_nb_epochs;
    double rate = (argc > 6) ? std::atof(argv[6]) : default_learning_rate;
    bool print = (argc > 7) ? std::atoi(argv[7]) == 1 : false;

    cerr << "If you see an assert failure now that probably means that" << endl
         << "\tsome of the neurons are not properly initialised or not properly connected."
         << endl;
    OneLayerPerceptron *perceptron = new OneLayerPerceptron(dim, size, rate, rate / epochs, sigma, sigma_der);

    cout << "Initialised a 1-layer perceptron" << endl;
    cout << "\tSize of the hidden layer:\t"
         << perceptron->getNbNeurons() << endl;
    cout << "\tLearning rate:\t\t\t"
         << perceptron->getLearningRate() << endl;
    cout << "\tDecay:\t\t\t\t" << perceptron->getDecay() << endl;

    cout << "Training the perceptron over "
         << epochs << " epochs for regression over column "
         << regr << endl;

    run_on_data(perceptron, training, epochs, regr, true, print);

    cout << "Switching off learning...\t";

    perceptron->initLearningRate(0);

    cout << "done. Learning rate = "
         << perceptron->getLearningRate() << endl;

    cout << "Testing the perceptron on the training data "
         << "(" << epochs << " times)" << endl;

    run_on_data(perceptron, training, epochs, regr, false, print);

    cout << "Testing the perceptron on the testing data "
         << "(" << epochs << " times)" << endl;

    run_on_data(perceptron, testing, epochs, regr, false, print);

    cout << "Deleting the perceptron...\t";
    delete perceptron;
    cout << "done." << endl;
    return 0;
}

int main(int argc, char *argv[])
{
    if (argc > 1)
    {
        if (argc < 3)
        {
            cout << "Usage: " << endl
                 << argv[0] << " <train_file> <test_file> " << endl
                 << "\t[ <column_for_regression> [ <nuber_of_neurons> [ <nuber_of_epochs> [ <learning_rate> [ <print> ] ] ] ] ] " << endl;
            cout << argv[0] << endl
                 << endl;

            cout << "The second form is for unitary testing" << endl;
            return 1;
        }

        return run(argc, argv);
    }
    else 
        return (
            test3() 
            && 
            test4("../csv/test.csv")
        ) ? 0 : 1;
}
