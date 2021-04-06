
#include "point.hpp"
#include "cloud.hpp"
#include "dendrogram.hpp"

#include <cstdlib> // for rand, srand
#include <ctime>   // for time
#include <cassert>
#include <cmath>   // for abs

#include <iostream>
using std::cout;
using std::endl;

void print_point(dendrogram *d, int i) {
    cout << i << " (" << d->get_name(i) << ")";
}

int trace_find(dendrogram *d, int i) {
    int p = d->get_parent(i);

    cout << d->get_rank(i) << ", ";
    print_point(d, i);
    cout << ", ";
    if (p != -1) {
        print_point(d, p);
        cout << ", " << d->get_height(i);
    }
    cout << endl;

    if (p == -1)
        return i;

    return trace_find(d, p);
}

void init_cloud(cloud &c, double *points, std::string *labels, int size) {
    cout << "with points ";
    // temporary container
    point p;
    for (int i = 0; i < size; i++) {
        cout << points[i] << " ";
        p.coords[0] = points[i];
        p.name = labels[i];
        c.add_point(p);
    }
    cout << "\t";
}

class test_dendro : public dendrogram {
public:
    // Memory leak here, yeah, to fix
    test_dendro (cloud &_c) : dendrogram (*(new graph(_c))) {}

    void merge (int p1, int p2, double d) {
        dendrogram::merge(new edge (p1, p2, d));
    }

    bool test_find (int *truth) {
        for (int i = 0; i < g->get_num_nodes(); i++) {
            if (find(i) != truth[i])
                return false;
        }
        return true;
    }

    bool test_parent (int *truth) {
        for (int i = 0; i < g->get_num_nodes(); i++) {
            if (parent[i] != truth[i])
                return false;
        }
        return true;
    }

    bool test_rank(int *truth) {
        for (int i = 0; i < g->get_num_nodes(); i++) {
            if (rank[i] != truth[i])
                return false;
        }
        return true;
    }

    bool test_left(int *truth) {
        for (int i = 0; i < g->get_num_nodes(); i++) {
            if (left[i] != truth[i])
                return false;
        }
        return true;
    }

    bool test_down(int *truth) {
        for (int i = 0; i < g->get_num_nodes(); i++) {
            if (down[i] != truth[i])
                return false;
        }
        return true;
    }

    bool test_height(double *truth) {
        for (int i = 0; i < g->get_num_nodes(); i++) {
            if (height[i] != truth[i])
                return false;
        }
        return true;
    }

    bool test_clusters(int *truth) {
        for (int i = 0; i < g->get_num_nodes(); i++) {
            if (clusters[i] != truth[i])
                return false;
        }
        return true;
    }

    bool test_get_cluster_height (double *truth, int *true_clusters) {
        for (int i = 0; i < g->get_num_nodes(); i++) 
            if (i == true_clusters[i] && std::abs(get_cluster_height(i) - truth[i]) > 0.01)
                    return false;                    
        return true;
    }

    graph *get_graph() {
        return g;
    }
};

std::string *test_find () {
    cloud c (1, 6);
    double points[] = {0.0, 1.0, 3.0, 4.0, 9.0};
    std::string labels[] = {"0", "1", "3", "4", "9"};
    init_cloud(c, points, labels, 5);

    test_dendro dendro (c);
    int find[] = {0, 1, 2, 3, 4};
    int parent[] = {-1, -1, -1, -1, -1};
    if (!dendro.test_find(find))
        return new std::string ("incorrect find() on initial dendrogram");

    parent[0] = 1;
    dendro.set_parent(parent);
    find[0] = 1;
    if (!dendro.test_find(find))
        return new std::string ("incorrect find() after one-level merge");

    parent[2] = 3;
    dendro.set_parent(parent);
    find[2] = 3;
    if (!dendro.test_find(find))
        return new std::string("incorrect find() after one-level merge");

    parent[1] = 2;
    dendro.set_parent(parent);
    find[0] = 3;
    find[1] = 3;
    if (!dendro.test_find(find))
        return new std::string("incorrect find() after two-level merge");;

    parent[4] = 3;
    dendro.set_parent(parent);
    find[4] = 3;
    if (!dendro.test_find(find))
        return new std::string("incorrect find() on final test dendrogram");

    return new std::string("[OK]");
}

std::string *test_merge () {
    cloud c (1, 5);
    double points[] = {0.0, 1.0, 3.0, 4.0, 9.0};
    std::string labels[] = {"0", "1", "3", "4", "9"};
    init_cloud(c, points, labels, 5);

    test_dendro dendro (c);
    dendro.merge(1, 0, 1);
    dendro.merge(3, 2, 1);
    dendro.merge(2, 1, 2);
    dendro.merge(4, 3, 5);

    int parent[] = {  1,  3,   3,  -1,   3};
    int rank[]   = {  0,  1,   0,   2,   0};
    int left[]   = { -1,  2,  -1,  -1,   1};
    int down[]   = { -1,  0,  -1,   4,  -1};
    double height[] = {0.5,  1, 0.5,  -1, 2.5};

    if (!dendro.test_parent(parent))
        return new std::string ("incorrect parent after merge");

    if (!dendro.test_rank(rank))
        return new std::string ("incorrect rank after merge");

    if (!dendro.test_left(left))
        return new std::string ("incorrect left after merge");

    if (!dendro.test_down(down))
        return new std::string ("incorrect rank after merge");

    if (!dendro.test_height(height))
        return new std::string ("incorrect height after merge");

    return new std::string ("[OK]");
}

std::string *test_build () {
    cloud c (1, 5);
    double points[] = {0.0, 1.0, 3.0, 4.0, 9.0};
    std::string labels[] = {"0", "1", "3", "4", "9"};
    init_cloud(c, points, labels, 5);

    test_dendro dendro (c);
    dendro.build();

    int parent[] = {  1,  3,   3,  -1,   3};
    int rank[]   = {  0,  1,   0,   2,   0};
    int left[]   = { -1,  2,  -1,  -1,   1};
    int down[]   = { -1,  0,  -1,   4,  -1};
    double height[] = {0.5,  1, 0.5,  -1, 2.5};

    int parent_alt[] = { -1,  0,  0,  2,   0};
    int rank_alt[]   = {  2,  0,  1,   0,  0};
    int left_alt[]   = { -1, -1,  1,  -1,   2};
    int down_alt[]   = {  4, -1,  3,  -1,  -1};
    double height_alt[] = {-1, 0.5, 1, 0.5, 2.5};

    if (!dendro.test_parent(parent) && !dendro.test_parent(parent_alt))
        return new std::string ("incorrect parent after build");

    if (!dendro.test_rank(rank) && !dendro.test_rank(rank_alt))
        return new std::string ("incorrect rank after build");

    if (!dendro.test_left(left) && !dendro.test_left(left_alt))
        return new std::string ("incorrect left after build");

    if (!dendro.test_down(down) && !dendro.test_down(down_alt))
        return new std::string ("incorrect rank after build");

    if (!dendro.test_height(height) && !dendro.test_height(height_alt))
        return new std::string ("incorrect height after build");

    return new std::string ("[OK]");
}

std::string *test_set_clusters() {
    cloud c (1, 5);
    double points[] = {0.0, 1.0, 3.0, 4.0, 9.0};
    std::string labels[] = {"0", "1", "3", "4", "9"};
    init_cloud(c, points, labels, 5);

    test_dendro dendro (c);
    dendro.build();

    int clusters[] = {-1, -1, -1, -1, -1};
    int clusters_alt[] = {-1, -1, -1, -1, -1};
    if (!dendro.test_clusters(clusters))
        return new std::string ("incorrect clusters upon initialisation");

    dendro.clear_clusters();
    dendro.set_clusters(0.5);
    clusters[0] = 1;
    clusters[1] = 1;
    clusters[2] = 3;
    clusters[3] = 3;
    clusters[4] = 4;
    if (!dendro.test_clusters(clusters))
        return new std::string ("incorrect clusters at height 0.5");

    dendro.clear_clusters();
    dendro.set_clusters(0.7);
    if (!dendro.test_clusters(clusters))
        return new std::string ("incorrect clusters at height 0.7");

    dendro.clear_clusters();
    dendro.set_clusters(1.0);
    clusters[0] = 3;
    clusters[1] = 3;
    clusters[2] = 3;
    clusters[3] = 3;
    clusters[4] = 4;
    if (!dendro.test_clusters(clusters))
        return new std::string ("incorrect clusters at height 1");

    dendro.clear_clusters();
    dendro.set_clusters(2.5);
    clusters[4] = 3;
    if (!dendro.test_clusters(clusters))
        return new std::string ("incorrect clusters at height 2.5");

    return new std::string ("[OK]");
}

std::string *test_count_clusters () {
    cloud c (1, 5);
    double points[] = {0.0, 1.0, 3.0, 4.0, 9.0};
    std::string labels[] = {"0", "1", "3", "4", "9"};
    init_cloud(c, points, labels, 5);

    test_dendro dendro (c);
    dendro.build();

    int clusters[] = {-1, -1, -1, -1, -1};
    if (!dendro.test_clusters(clusters))
        return new std::string ("incorrect clusters upon initialisation");

    dendro.clear_clusters();
    dendro.set_clusters(0.5);
    if (dendro.count_ns_clusters() != 2) {
        return new std::string ("incorrect number of clusters at height 0.5 -- should be 2");
    }

    dendro.clear_clusters();
    dendro.set_clusters(0.7);
    if (dendro.count_ns_clusters() != 2)
        return new std::string ("incorrect number of clusters at height 0.7 -- should be 2");

    dendro.clear_clusters();
    dendro.set_clusters(1.0);
    if (dendro.count_ns_clusters() != 1)
        return new std::string ("incorrect number of clusters at height 1 -- should be 1");

    dendro.clear_clusters();
    dendro.set_clusters(2.5);
    if (dendro.count_ns_clusters() != 1)
        return new std::string ("incorrect number of clusters at height 2.5 -- should be 1");

    return new std::string ("[OK]");
}

std::string *test_get_cluster_height () {
    cloud c (1, 5);
    double points[] = {0.0, 1.0, 3.0, 4.0, 9.0};
    std::string labels[] = {"0", "1", "3", "4", "9"};
    init_cloud(c, points, labels, 5);

    test_dendro dendro (c);
    dendro.build();

    dendro.clear_clusters();
    dendro.set_clusters(0.5);
    double height[] = {0,  0.5, 0,  0.5, 0};
    int clusters[] = {1, 1, 3, 3, 4};
    double height_alt[] = {0.5, 0, 0.5, 0, 0};
    int clusters_alt[] = {0, 0, 2, 2, 4};
    if (!dendro.test_get_cluster_height(height, clusters) && !dendro.test_get_cluster_height(height_alt, clusters_alt))
        return new std::string ("cluster heights incorrect after setting clusters at height 0.5");

    dendro.clear_clusters();
    dendro.set_clusters(1);
    height[1] = 0;
    height[3] = 1;
    clusters[0] = 3;
    clusters[1] = 3;
    height_alt[0] = 1;
    height_alt[2] = 0;
    clusters_alt[2] = 0;
    clusters_alt[3] = 0;
    if (!dendro.test_get_cluster_height(height, clusters) && !dendro.test_get_cluster_height(height_alt, clusters_alt))
        return new std::string ("cluster heights incorrect after setting clusters at height 1");

    dendro.clear_clusters();
    dendro.set_clusters(2.5);
    height[3] = 2.5;
    clusters[4] = 3;
    height_alt[0] = 2.5;
    clusters_alt[4] = 0;
    if (!dendro.test_get_cluster_height(height, clusters) && !dendro.test_get_cluster_height(height_alt, clusters_alt))
        return new std::string ("cluster heights incorrect after setting clusters at height 2.5");

    return new std::string ("[OK]");
}

int main(int argc, char *argv[]) {
    srand(time(NULL));

    if (argc < 2) {
        cout << "Usage (clouds): " << endl
             << argv[0] << " " << "<file name> <data dimension> <nb points>" << endl
             << "Example: " << argv[0] << " ../../data/iris.data 4 150" << endl << endl
             << "Usage (distance matrix):" << argv[0] << " "
             << "<file name>" << endl
             << "Example: " << argv[0] << " ../../data/test5_matrix.data" << endl << endl
             << "Usage (testing): " << endl
             << argv[0] << " test" << endl;
        return 0;
    }

    std::ifstream is(argv[1]);
    if (is.is_open()) {
        graph* g;
        if (argc == 4) {
            int d = std::stoi(argv[2]);
            int nmax = std::stoi(argv[3]);

            cloud c(d, nmax);

            c.load(is);

            cout << "Loaded "
                << c.get_n()
                << " points from "
                << argv[1];
            cout << ((c.get_n() == nmax) ? "\t[OK]" : "\t[NOK]") << endl;
            if (c.get_n() != nmax)
                cout << "Must have loaded " << nmax << " points" << endl;
            g = new graph(c);
        } else {
            g = graph::load_matrix(is);
            cout << "Loaded "
                 << g->get_num_nodes()
                 << " vertices from "
                 << argv[1] << endl;    
        }
        is.close();
    
        dendrogram dg(*g);
        dg.build();

        cout << "Height of the dendrogram:\t"
            << dg.get_dendro_height() << endl
            << "(For iris.data, height should be 0.820061)" << endl;

        cout << "Printing traces to root from 10 random points..." << endl;
        cout << "Rank, point, parent, height" << endl;
        srand(static_cast<unsigned>(time(0)));
        for (int i = 0; i < 10; i++) {
            int node = rand() % dg.get_n();
            trace_find(&dg, node);
            cout << endl;
        }

        double eps = 0.01;
        cout << "Looking for significant heights (up to " << eps << ")...";
        dg.find_heights(eps);
        cout << "\tdone" << endl;
        int count = dg.get_count_sign_heights();
        cout << "Number of significant heights (up to "
            << eps << ") = " << count << ": " << endl;
        for (int i = 0; i < count; i++)
            cout << dg.get_sign_height(i) << " ";
        cout << endl;

        cout << "Printing clusters at significant heights" << endl;
        for (int i = 0; i < count; i++) {
            dg.clear_clusters();
            dg.set_clusters(dg.get_sign_height(i));
            cout << dg.count_ns_clusters()
                << " cluster"
                << (dg.count_ns_clusters() > 1 ? "s" : "")
                << " found at height " << dg.get_sign_height(i) << endl;
            dg.print_clusters();
            dg.print_dendrogram();
            cout << endl;
        }
    } else {
        std::string *msg;

        cout << "Testing dendrogramm::find()...\t";
        msg = test_find();
        cout << *msg << endl;
        delete msg;

        cout << "Testing dendrogramm::merge()...\t";
        msg = test_merge();
        cout << *msg << endl;
        delete msg;

        cout << "Testing dendrogramm::build()...\t";
        msg = test_build();
        cout << *msg << endl;
        delete msg;        

        cout << "Testing dendrogram::set_clusters()...\t";
        msg = test_set_clusters();
        cout << *msg << endl;
        delete msg;

        cout << "Testing dendrogram::count_clusters()...\t";
        msg = test_count_clusters();
        cout << *msg << endl;
        delete msg;

        cout << "Testing dendrogram::get_cluster_height()...\t";
        msg = test_get_cluster_height();
        cout << *msg << endl;
        delete msg;        
    }

    return 0;
}
