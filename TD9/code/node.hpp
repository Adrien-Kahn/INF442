#pragma once

#include <iostream>

class Node
{
private:
    // Forward-propagated signal
    double signal;

    // Back-propagated error
    // Must be initialised to 0 (see comment to step_back() )
    double back_value;

    // For testing only
    static int count;
public:
    Node();
    Node(double _signal);
    ~Node();
    
    double getSignal() const;
    double getBackValue() const;

    void setSignal(double _signal);
    void setBackValue(double _back_value);

    static int getCount();
};

std::ostream &operator<<(std::ostream &str, const Node &n);
