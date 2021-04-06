#pragma once

#include "cloud.hpp"

/* edge -- pairs of cloud point indices with methods to
 *  - compute the distance from source to target
 *  - compare the lengths of two edges -- needed for sorting
 */

class edge {
    int p1, p2;
    double length;

    static int count_compare; // For testing only

public:
    edge(int _p1, int _p2, double _length);
    ~edge();
    
    int get_p1();
    int get_p2();
    double get_length();

    static bool compare(edge *e1, edge *e2);
    static int get_count_compare();
};
