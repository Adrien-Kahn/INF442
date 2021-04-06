#pragma once

#include "cloud.hpp"
#include "graph.hpp"

class dendrogram {
protected:
    graph *g;

    // Union-find structure elements (exercise 3.2)
    int *parent;
    int *rank;

    // Binary tree structure elements (exercise 3.3)
    int *left;
    int *down;

    // Height in the dendrogram at which
    // the represented cluster is merged into the parent one
    double *height;

    // Clusters
    double cut_height;
    int *clusters; // List of lists: each cluster is a list with representatives being heads
    int total_clusters; // Number of clusters currently identified
    int nonsingle_clusters; // Number of nonsingleton clusters

    // List of significant heights (exercise 3.4)
    double *sign_heights = nullptr;
    int count_sign_heights;

    void merge(edge *e);
    int _count_ns_clusters(); // counts the nonsingular clusters from scratch

public:
    dendrogram(graph& _g);
    ~dendrogram();
    void build();

    int find(int i);            // Find the cluster representative in the union-find structure
    double get_dendro_height(); // Get the height of the dendrogram

    void find_heights(double eps); // Put all heights up to eps into sign_heights

    void set_clusters(int i, double h);
    void set_clusters(double h);
    int count_ns_clusters();
    void clear_clusters();

    double get_cluster_height(int cluster);

    /*** GETTERS ***/
    int get_parent(int i);
    int get_rank(int i);
    int get_cluster(int i);
    double get_height(int i);
    std::string& get_name(int i);

    int get_down(int i);
    int get_left(int i);

    int get_count_sign_heights();
    double get_sign_height(int i);
    double get_cut_height();
    int get_total_clusters();
    long get_n();

    // For testing only
    void print_node(int i);
    void print_dendrogram();
    void print_clusters();
    void iterate_sign_heights();
    void set_parent(int *_parent);
};
