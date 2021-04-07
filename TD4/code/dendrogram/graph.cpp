#include "graph.hpp" // This is the header for the class implemented here

#include "cloud.hpp" // Used in the constructor
#include "edge.hpp"  // Used in almost all methods

#include <algorithm> // This provides the sort() method for the constructor
#include <cassert>
#include <iostream>

/* graph -- method implementations */

graph::graph(cloud &_c) {
    n = _c.get_n();
    size = n * (n - 1) / 2;
    
    edges = new edge* [size];
    edge **edge_iterator = edges;
    
    for (int i = 0; i < n; i++) {
    	for (int j = i + 1; j < n; j++) {
    		edge *e = new edge(i, j, _c.get_point(i).dist(_c.get_point(j)));
    		*edge_iterator = e;
    		edge_iterator++;
    	}
    }
    
    std::sort(edges, edges + size, edge::compare);
    
    iterator_pos = 0;
}

graph::graph(long _n, std::string _node_names[], double** dist_matrix) {
    n = _n;
    size = n * (n - 1) / 2;
    node_names = new std::string[n];
    node_names = _node_names;
    
    edges = new edge* [size];
    edge **edge_iterator = edges;
    
    for (int i = 0; i < n; i++) {
    	for (int j = i + 1; j < n; j++) {
    		edge *e = new edge(i, j, dist_matrix[i][j]);
    		*edge_iterator = e;
    		edge_iterator++;
    	}
    }
    
    std::sort(edges, edges + size, edge::compare);


    iterator_pos = 0;

}

graph::~graph() {
    
    edge *e;
    start_iteration();
    while ((e = get_next()) != NULL) {
        delete e;
    }
    delete[] edges;
    //delete[] node_names;
}

long graph::get_num_edges() {
    return size;
}

std::string& graph::get_name(int i) {
    assert(i >= 0 && i < n);
    return node_names[i];
}

long graph::get_num_nodes() {
    return n;
}

void graph::start_iteration() {
    iterator_pos = 0;
}

edge *graph::get_next() {
    
    if (iterator_pos == size) return NULL;
    
    edge *eptr = *(edges + iterator_pos);
    iterator_pos++;
    
    return eptr;
}

graph* graph::load_matrix(std::ifstream &is) {
    assert(is.is_open());

    int n;
    is >> n;
    std::string node_names[n];
    for (size_t i = 0; i < n; ++i) {
        is >> node_names[i];
    }

    std::cout << "Names read" << std::endl;
    double** dist_matrix = new double*[n];
    for (size_t i = 0; i < n; ++i) {
        dist_matrix[i] = new double[n];
        for (size_t j = 0; j < n; ++j) {
            is >> dist_matrix[i][j];
        }
    }

    graph* g = new graph(n, node_names, dist_matrix);
    for (size_t i = 0; i < n; ++i) {
        delete[] dist_matrix[i];
    }
    delete[] dist_matrix;

    return g;
}
