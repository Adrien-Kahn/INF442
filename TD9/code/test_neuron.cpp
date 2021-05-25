#include "neuron.hpp"

#include <iostream>
#include <memory>
#include <random>
#include <mutex>
#include <atomic>

int main()
{
    using std::cout;
    using std::endl;

    cout << "Testing Neuron::step()...\t";
    std::vector<Node *> inputs({new Node(), new Node(), new Node()});
    Neuron neuron(
        3, inputs,
        [](double x) -> double { return x; },
        [](double x) -> double { return 1; });
    neuron.setWeight(0, 1);
    neuron.setWeight(1, 10);
    neuron.setWeight(2, 100);
    neuron.setBias(1);

    Neuron neuron2(
        3, inputs,
        [](double x) -> double { return x + 1; },
        [](double x) -> double { return 1; });
    neuron2.setWeight(0, 1);
    neuron2.setWeight(1, 10);
    neuron2.setWeight(2, 100);
    neuron2.setBias(1);

    inputs[0]->setSignal(10);
    inputs[1]->setSignal(100);
    inputs[2]->setSignal(1000);

    // cout << "At neuron initialisation" << endl;
    // cout << neuron << endl;
    // cout << neuron2 << endl;

    neuron.step();
    neuron2.step();
    // cout << "After one step" << endl;
    // cout << neuron << endl;
    // cout << neuron2 << endl;
    int n1s1 = neuron.getAxon()->getSignal();
    int n2s1 = neuron2.getAxon()->getSignal();
    int ci1s1 = neuron.getCollectedInput();
    int ci2s1 = neuron2.getCollectedInput();

    neuron.step();
    neuron2.step();
    // cout << "After two steps" << endl;
    // cout << neuron << endl;
    // cout << neuron2 << endl;
    int n1s2 = neuron.getAxon()->getSignal();
    int n2s2 = neuron2.getAxon()->getSignal();
    int ci1s2 = neuron.getCollectedInput();
    int ci2s2 = neuron2.getCollectedInput();

    bool ex1ok = false;

    if (n2s1 != n1s1 + 1 || n2s2 != n1s2 + 1)
        cout << "[NOK] - Activation function not used correctly" << endl;
    else  if (n1s1 != n1s2 || n2s1 != n2s2)
        cout << "[NOK] - Collected input not initialised properly" << endl;
    else if (ci1s1 != n1s1 || ci1s2 != n1s2)
        cout << "[NOK] - Collected input not computed properly (did you use a local variable instead of Neuron::collected_input?)" << endl;
    else if (ci2s1 == n2s1)
        cout << "[NOK] - Collected input not computed properly (did you apply the activation function to the stored value?)" << endl;
    else  switch (n1s1)
    {
    case 101009:
        cout << "[OK]" << endl;
        ex1ok = true;
        break;

    case 101010:
        cout << "[NOK] - Bias not used properly" << endl;
        break;

    case 1109:
        cout << "[NOK] - Weights not used properly" << endl;
        break;

    case 1110:
        cout << "[NOK] - Bias and weights are not used properly" << endl;
        break;

    case 1009:
        cout << "[NOK] - Last dendrite not used properly" << endl;
        break;

    case 112:
        cout << "[NOK] - Signals not used properly" << endl;
        break;

    case 111:
        cout << "[NOK] - Bias and signals not used properly" << endl;
        break;

    default:
        cout << "[NOK] - Probably more than one error" << endl;
        break;
    }

    if (!ex1ok)
    {
        delete inputs[0];
        delete inputs[1];
        delete inputs[2];
        return 0;
    }

    cout << "Testing Neuron::step_back()...\t";

    neuron.setWeight(0, 1);
    neuron.setWeight(1, 10);
    neuron.setWeight(2, 100);
    neuron.setBias(1);

    std::vector<Node *> inputs3({new Node(), new Node(), new Node()});
    inputs3[0]->setSignal(10);
    inputs3[1]->setSignal(100);
    inputs3[2]->setSignal(1000);

    Neuron neuron3(
        3, inputs3,
        [](double x) -> double { return x*x; },
        [](double x) -> double { return 2*x; });
    neuron3.setWeight(0, 1);
    neuron3.setWeight(1, 10);
    neuron3.setWeight(2, 100);
    neuron3.setBias(1);

    neuron.setCollectedInput(5);
    neuron3.setCollectedInput(6);
    neuron.setLearningRate(3);
    neuron3.setLearningRate(3);

    neuron.setBackValue(2);
    neuron3.setBackValue(2);

    // cout << "At neuron initialisation" << endl;
    // cout << neuron << endl;
    // cout << neuron3 << endl << endl;

    inputs[2]->setBackValue(-2);
    neuron.step_back();
    neuron3.step_back();

    // cout << "After one step backwards" << endl;
    // cout << neuron << endl;
    // cout << neuron3 << endl << endl;

    int n1w1[4];
    int n1back1[4];
    int n3w1[4];
    int n3back1[4];
    for (int i = 0; i < 3; i++) {
        n1w1[i+1] = neuron.getWeight(i);
        n3w1[i+1] = neuron3.getWeight(i);
        n1back1[i+1] = neuron.getDendrite(i)->getBackValue();        
        n3back1[i+1] = neuron3.getDendrite(i)->getBackValue();        
    }
    n1w1[0] = neuron.getBias();
    n3w1[0] = neuron3.getBias();
    n1back1[0] = neuron.getBiasDendrite()->getBackValue();        
    n3back1[0] = neuron3.getBiasDendrite()->getBackValue();        

    if (n1w1[3] == 100) 
        cout << "[NOK] - Last dendrite not updated properly" << endl;
    else if (n1w1[0] == 1)
        cout << "[NOK] - Bias not updated properly" << endl;
    else if (n1back1[1] == 0)
        cout << "[NOK] - Error not propagated to next level" << endl;
    else if (n1back1[2] == 2 || n3back1[2] == 24)
        cout << "[NOK] - Forgot the weights during error propagation?" << endl;
    else if (n1back1[2] == 10 || n1back1[2] == 100)
        cout << "[NOK] - Forgot the axon back value during error propagation?" << endl;
    else if (n1back1[3] == 198)
        cout << "[NOK] - Do you add propagated error to the existing one?" << endl;
    else if (n1back1[3] == -202)
        cout << "[NOK] - Do you subtract propagated error from the existing one?" << endl;
    else if (n1w1[1] == 60 && n1w1[2] == 600 && n1w1[3] == 6000)
        cout << "[NOK] - '=' instead of '-=' in weight update?" << endl;
    else if (n1w1[1] == -60 && n1w1[2] == -600 && n1w1[3] == -6000)
        cout << "[NOK] - '= -' (assignment) instead of '-=' (cumulative updtate) in weight update?" << endl;
    else if (n1w1[1] == -19 && n1w1[2] == -190 && n1w1[3] == -1900)
        cout << "[NOK] - Forgot the learning rate in weight update?" << endl;
    else if (n1w1[1] == -29 && n1w1[2] == -290 && n1w1[3] == -2900)
        cout << "[NOK] - Forgot the axon back value in weight update?" << endl;
    else if (n1w1[1] == -5 && n1w1[2] == 4 && n1w1[3] == 94)
        cout << "[NOK] - Forgot the dendrite signal in weight update?" << endl;
    else if (n3w1[1] == -59 && n3w1[2] == -590 && n3w1[3] == -5900)
        cout << "[NOK] - Forgot the collected input in weight update?" << endl;
    else if (n3w1[1] == -359 && n3w1[2] == -3590 && n3w1[3] == -35900)
        cout << "[NOK] - Forgot the activation function derivative in weight update?" << endl;
    else if (n3w1[1] == -2159 && n3w1[2] == -21590 && n3w1[3] == -215900)
        cout << "[NOK] - Used the activation function instead of its derivative in weight update?" << endl;
    else if ((n1w1[0] == 7 && n1w1[1] == -59 && n1w1[2] == -590 && n1w1[3] == -5900) &&
             (n3w1[0] == 73 && n3w1[1] == -719 && n3w1[2] == -7190 && n3w1[3] == -71900) &&
             ((n1back1[0] == 2 || n1back1[0] == 0) && n1back1[1] == 2 && n1back1[2] == 20 && n1back1[3] == 200) &&
             ((n3back1[0] == 24 || n3back1[0] == 0) && n3back1[1] == 24 && n3back1[2] == 240 && n3back1[3] == 2400))
        cout << "[OK]" << endl;
    else {
        cout << "[NOK] - Probably more than one error" << endl;
    }

    delete inputs[0];
    delete inputs[1];
    delete inputs[2];
    delete inputs3[0];
    delete inputs3[1];
    delete inputs3[2];

    return 0;
}
