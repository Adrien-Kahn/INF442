#ifndef RETRIEVAL_HPP
#define RETRIEVAL_HPP
#include <algorithm>  // for sort
#include <cassert>    // for assertions
#include <cfloat>     // for DBL_MAX
#include <cmath>      // for math operations like sqrt, log, etc
#include <cstdlib>    // for rand, srand
#include <ctime>      // for clock
#include <fstream>    // for ifstream
#include <iostream>   // for cout

using std::cout;
using std::endl;

/*****************************************************
 * TD 2: K-Dimensional Tree (kd-tree)                *
 *****************************************************/

const bool debug = false;  // debug flag, r

typedef double *point;  // point = array of coordinates (doubles)

void print_point(point p, int dim) {
    std::cout << "[ ";
    for (int i = 0; i < dim; ++i) {
        std::cout << p[i] << " ";
    }
    std::cout << "]" << std::endl;
}

void pure_print(point p, int dim) {
    std::cout << p[0];
    for (int j = 1; j < dim; j++)
        std::cout << " " << p[j];
    std::cout << "\n";
}


/**
 * This function computes the Euclidean distance between two points
 *
 * @param p  the first point
 * @param q the second point
 * @param dim the dimension of the space where the points live
 * @return the distance between p and q, i.e., the length of the segment pq
 */
double dist(point p, point q, int dim) {
    // Exercise 1.1
    double distance = 0.0;

    for (int k = 0; k < dim; k++) {
    	distance += (p[k] - q[k])*(p[k] - q[k]);
    }
    
    distance = sqrt(distance);

    return distance;
}

/**
 * This function for a given query point q  returns the index of the point
 * in an array of points P that is closest to q using the linear scan algorithm.
 *
 * @param q the query point
 * @param dim the dimension of the space where the points live
 * @param P a set of points
 * @param n the number of points in P
 * @return the index of the point in p that is closest to q
 */
int linear_scan(point q, int dim, point *P, int n) {
    // Exercise 1.2
    int idx = -1;  // It will contain the index of the closest point

    double mindist = dist(P[0], q, dim) + 1;

    for (int k = 0; k < n; k++) {
    	if (dist(P[k], q, dim) < mindist) {
    		mindist = dist(P[k], q, dim);
    		idx = k;
    	}
    }

    return idx;
}

/**
 * This function computes the median of all the c coordinates 
 * of an subarray P of n points that is P[start] .. P[end - 1]
 *
 * @param P a set of points
 * @param start the starting index
 * @param end the last index; the element P[end] is not considered
 * @return the median of the c coordinate
 */

// Auxiliary function gets the index of the minimum of P along the axis c within the start-end range (end excluded)
int minidx(point *P, int start, int end, int c) {
	int idx = start;
	for (int j = start; j < end; j++) {
		if (P[j][c] < P[idx][c]) {
		idx = j;
		}
	}
	return idx;
}


// Auxiliary function that selection sorts P between start and end (excluded) along the axis c in place
void sort(point *P, int start, int end, int c) {
	for (int i = start; i < end; i++) {
		int j = minidx(P, i, end, c);
		std::swap(P[i],P[j]);
	}
}

// Prints whether the array of points P is sorted along the axis c from start to end
void isSorted(point *P, int start, int end, int c) {
	for (int k = start; k < end - 1; k++) {
		if (P[k][c] > P[k + 1][c]) {
			std::cout << "Not sorted" << std::endl;
			return;
		}
	}
	std::cout << "Sorted" << std::endl;
	return;
}


double compute_median(point *P, int start, int end, int c) {
    sort(P, start, end, c);
    isSorted(P, start, end, c);
    return P[start + ((end - start)/2)][c];
}

double compute_median_1(point *P, int start, int end, int c) {
    double a[end - start];
    std::cout << "Hello there" << std::endl;
    print_point(P[end - 1], c);
    std::cout << "start: " << start << "    end: " << end << std::endl;
    if (end == 150) {
    	std::cout << "hallo" << std::endl;
    	print_point(P[149], 1);
    	std::cout << "end hallo" << std::endl;
    }
    int n = sizeof(a)/sizeof(a[0]);
    for (int k = start; k < end; k++) {
//    	std::cout << k << std::endl;
    	a[k] = P[k][c];
//    	std::cout << k << std::endl;
    }
    std::sort(a, a + n);
    return a[start + ((end - start)/2)];
}

/**
 * Partitions the the array P wrt to its median value of a coordinate
 *
 * @param P a set of points (an array)
 * @param start the starting index
 * @param end the last index; the element P[end] is not considered
 * @param c the coordinate that we will consider the median
 * @param dim the dimension where the points live
 * @return the index of the median value
 */
int partition(point *P, int start, int end, int c, int dim) {
    // Exercise 3.2
    double m = compute_median(P, start, end, c);
    int idx = -1;  // this is where we store the index of the median

    // TODO

    return idx;
}

typedef struct node {  // node for the kd-tree
    // Exercise 3.3

    // TODO
    int c;               // coordinate for split
    double m;            // split value
    int idx;             // index of data point stored at the node
    node *left, *right;  // children

} node;

/**
 * Creates a leaf node in the kd-tree
 *
 * @param val the value of the leaf node
 * @return a leaf node that contains val
 */
node *create_node(int _idx) {
    // Exercise 3.3

    // TODO
    
    // Do not forget to replace this return by a correct one!
    return new node;
}

/**
 * Creates a internal node in the kd-tree
 *
 * @param idx the value of the leaf node
 * @return an internal node in the kd-tree that contains val
 */
node *create_node(int _c, double _m, int _idx,
                  node *_left, node *_right) {  
    // Exercise 3.3

    // TODO

    // Do not forget to replace this return by a correct one!
    return new node;
}

node* build(point *P, int start, int end, int c, int dim) {
    // builds tree for sub-cloud P[start -> end-1]
    assert(end - start >= 0);
    if (debug)
        std::cout << "start=" << start << ", end=" << end << ", c=" << c
                  << std::endl;
    if (end - start == 0)  // no data point left to process
        return NULL;
    else if (end - start == 1)  // leaf node
        return create_node(start);
    else {  // internal node
        if (debug) {
            std::cout << "array:\n";
            for (int i = start; i < end; i++)
                print_point(P[i], dim);
            // std::cout << P[i] << ((i==end-1)?"\n":" ");
        }
        // compute partition
        // rearrange subarray (less-than-median first, more-than-median last)
        int p = partition(P, start, end, c, dim);
        if (p == -1) { return NULL; }
        double m = P[p][c];
        // prepare for recursive calls
        int cc = (c + 1) % dim;  // next coordinate
        return create_node(c, m, p, build(P, start, p, cc, dim),
                           build(P, p + 1, end, cc, dim));
    }
}

/**
 *  Defeatist search in a kd-tree
 *
 * @param n the roots of the kd-tree
 * @param q the query point
 * @param dim the dimension of the points
 * @param P a pointer to an array of points
 * @param res the distance of q to its NN in P
 * @param nnp the index of the NN of q in P
 * @return a leaf node that contains val
 */
void defeatist_search(node *n, point q, int dim, point *P, double &res, int &nnp) {
    // Exercise 3.4
    // TODO
}

/**
 *  Backtracking search in a kd-tree
 *
 * @param n the roots of the kd-tree
 * @param q the query point
 * @param dim the dimension of the points
 * @param P a pointer to an array of points
 * @param res the distance of q to its NN in P
 * @param nnp the index of the NN of q in P
 * @return a leaf node that contains val
 */
void backtracking_search(node *n, point q, int dim, point *P, double &res, int &nnp) {
    // Exercise 3.5
    // TODO
}

#endif  // RETRIEVAL_HPP
