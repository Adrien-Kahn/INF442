#include "node.hpp"

#include <cassert>
#include <iostream>

Node::Node() : signal(0), back_value(0) { 
    count++; 
}
Node::Node(double _signal) : signal(_signal), back_value(0) { 
    count++; 
}
Node::~Node() { 
    count--; 
}

double Node::getSignal() const
{
    return signal;
}

double Node::getBackValue() const
{
    return back_value;
}

void Node::setSignal(double _signal)
{
    assert(this != (void *)0xffffffffffffffff);
    signal = _signal;
}

void Node::setBackValue(double _back_value)
{
    back_value = _back_value;
}

int Node::count = 0;

int Node::getCount() {
    return count;
}

std::ostream &operator<<(std::ostream &str, const Node &n)
{
    return str << "(" << n.getSignal() << ", " << n.getBackValue() << ")";
}
