#include <iostream>
#include <cassert>
#include <cfloat> // for DBL_MAX
#include <fstream>

#include "cloud.hpp"
#include "edge.hpp"
#include "graph.hpp"

using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cout << "Usage (clouds): " << argv[0] << " "
             << "<file name> <data dimension> <nb points>" << endl
             << "Example: " << argv[0] << " iris.data 4 150" << endl
             << "Usage (distance matrix):" << argv[0] << " "
             << "<file name>" << endl
             << "Example: " << argv[0] << " ../../data/test5_matrix.data" << endl;
        return 0;
    }

    const int d = (argc > 2) ? std::stoi(argv[2]) : 4;
    const int nmax = (argc > 3) ? std::stoi(argv[3]) : 150;

    graph* g;

    if (argc == 4) {
        const int d = std::stoi(argv[2]);
        const int nmax = std::stoi(argv[3]);

        cloud c(d, nmax);

        std::ifstream is(argv[1]);
        assert(is.is_open());

        c.load(is);
        is.close();

        cout << "Loaded "
             << c.get_n()
             << " points from "
             << argv[1];
        cout << ((c.get_n() == nmax) ? "\t[OK]" : "\t[NOK]") << endl;
        if (c.get_n() != nmax)
            cout << "Must have loaded " << nmax << " points" << endl;
        g = new graph(c);
    } else {
        std::ifstream is(argv[1]);
        assert(is.is_open());
        g = graph::load_matrix(is);
        is.close();
        cout << "Loaded "
             << g->get_num_nodes()
             << " vertices from "
             << argv[1] << endl;
    }
    
    cout << "Initialisation of the complete graph\t[OK]" << endl;
    cout << "Expected number of edges:\t" << g->get_num_edges() << endl;

    int count = 0;
    int count_compare = edge::get_count_compare();
    edge *e;

    g->start_iteration();
    while ((e = g->get_next()) != NULL) {
        count++;
    }
    cout << "Actual number of edges:\t" << count;
    assert(g->get_num_edges() == count);
    cout << "\t\t[OK]" << endl;
    cout << "Number of edge comparisons for sorting:\t"
         << count_compare << endl;

    g->start_iteration();
    edge *e1 = g->get_next();
    edge *e2;

    double min = e1->get_length();
    double max = min;

    g->start_iteration();
    while (e1 != NULL && (e2 = g->get_next()) != NULL) {
        if (edge::compare(e2, e1)) {
            cout << "Edge sorting failed!!" << endl;
            break;
        }

        double d = e2->get_length();
        if (d < min)
            min = d;

        if (d > max)
            max = d;

        e1 = e2;
    }

    if (e2 == NULL) {
        cout << "Edge sorting\t\t\t\t[OK]" << endl;
        cout << "Number of edge comparisons for testing:\t"
             << edge::get_count_compare() - count_compare << endl;

        cout << "Min distance between points:\t" << min << endl;
        cout << "Max distance between points:\t" << max << endl;
        cout << "Last distnace between points:\t" << e1->get_length() << endl;
    }

    delete g;

    return 0;
}
