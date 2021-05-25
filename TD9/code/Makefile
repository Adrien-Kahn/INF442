CXX = g++ 
CXXFLAGS = -g --std=c++11

.PHONY : all clean

all : test_neuron test_perceptron 

activation.o:	activation.cpp activation.hpp
dataset.o:		dataset.cpp dataset.hpp
node.o:			node.cpp node.hpp
neuron.o:		neuron.cpp neuron.hpp node.hpp
perceptron.o:	perceptron.cpp perceptron.hpp neuron.hpp dataset.hpp

test_neuron : test_neuron.cpp neuron.o node.o
	$(CXX) $(CXXFLAGS) -o $@ $^

test_perceptron : test_perceptron.cpp perceptron.o neuron.o node.o dataset.o activation.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean :
	rm -f test_neuron test_perceptron
	rm -f *~ *.o
	rm -rf *.dSYM
