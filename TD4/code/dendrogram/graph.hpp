#pragma once

#include "cloud.hpp"
#include "edge.hpp"

/*  graph -- an array of edge pointers arranged by
 *  increasing length.  Allows iteration:
 *
 *  init_iteration() places current at
 *    the start of the array
 *  get_next() returns the current edge
 *    and advances to the next one
 */
class graph {
    edge **edges;

    long n;
    long size;
    long iterator_pos;

    std::string* node_names;
public:
    graph(cloud &_c);
    graph(long _n, std::string _node_names[], double** dist_matrix);
    ~graph();

    void start_iteration();
    edge *get_next();

    long get_num_edges();
    long get_num_nodes();
    std::string& get_name(int i);

    static graph* load_matrix(std::ifstream &is);
};
