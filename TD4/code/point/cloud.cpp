#include "cloud.hpp"

#include <cassert>
#include <cfloat>	// for DBL_MAX
#include <cstdlib>	// for rand, srand


cloud::cloud(int _d, int _nmax, int _k) {
    d = _d;
    point::set_dim(_d);
    n = 0;
    k = _k;

    nmax = _nmax;

    points = new point[nmax];
    centers = new point[k];
}

cloud::~cloud() {
    delete[] centers;
    delete[] points;
}

void cloud::add_point(point &p, int label) {
    assert(n < nmax);

    for(int m = 0; m < d; m++) {
        points[n].coords[m] = p.coords[m];
    }
    points[n].label = label;

    n++;
}

int cloud::get_d() {
    return d;
}

int cloud::get_n() {
    return n;
}

int cloud::get_k() {
    return k;
}

point& cloud::get_point(int i) {
    return points[i];
}

point& cloud::get_center(int j) {
    return centers[j];
}

void cloud::set_center(point &p, int j) {
    for(int m = 0; m < d; m++)
        centers[j].coords[m] = p.coords[m];
}

double cloud::intracluster_variance() {
    double sum = 0.0;
    for(int i = 0; i < n; i++)
        sum += points[i].squared_dist(centers[points[i].label]);

    return sum / n;
}

int cloud::set_voronoi_labels() {
    int num_changes = 0;
    for (int i = 0; i < n; i++) {
        double min_dist = DBL_MAX;
        int min_ind = -1;

        for (int j = 0; j < k; j++) {
	        if (points[i].squared_dist(centers[j]) < min_dist) {
                min_dist = points[i].squared_dist(centers[j]);
                min_ind = j;
	        }
        }

        if (min_ind != points[i].label) {
            ++num_changes;
            points[i].label = min_ind;
        }
    }
    return num_changes;
}

void cloud::set_centroid_centers() {
    int *counts = new int[k];
    for (int j = 0; j < k; j++)
        counts[j] = 0;
    for (int i = 0; i < n; i++)
        counts[points[i].label]++;


    for (int j = 0; j < k; j++) {
        if(counts[j] != 0) { 
            for(int m = 0; m < d; m++) 
                centers[j].coords[m] = 0.0;
        }
    }

    for (int i = 0; i < n; i++) {
        for(int m = 0; m < d; m++)
            centers[points[i].label].coords[m] += points[i].coords[m];
    }

    for (int j = 0; j < k; j++) {
        if (counts[j] != 0) {
            for (int m = 0; m < d; m++)
                centers[j].coords[m] /= counts[j];
        }
    }

    delete[] counts;
}

void cloud::kmeans() {
    // set_centroid_centers(); // trivial initialization
    init_forgy();
    // init_plusplus();
    // init_random_partition();

    bool stabilized = false;

    while(true) {
        int num_changed = set_voronoi_labels();
        set_centroid_centers();
        if (num_changed == 0)
            break;
    }
}

void cloud::init_forgy() {
    int *already_chosen = new int[n];

    for (int j = 0; j < k; j++) {
        // choose index i different from those already chosen
        int i;
        bool new_index = false;
        while (!new_index) {
            i = rand() % n;
            // check whether i was already chosen
            new_index = true;
            for (int r = 0; r < j; r++) {
	            if(already_chosen[r] == i) {
                    new_index = false;
                    break;
                }
            }
        }
        already_chosen[j] = i;
        
        for (int m = 0; m < d; m++)
            centers[j].coords[m] = points[i].coords[m];
    }

    delete[] already_chosen;
}

void cloud::init_plusplus() {
    // choose first center
    int i = rand() % n;
    for (int m = 0; m < d; m++)
        centers[0].coords[m] = points[i].coords[m];

    // number of centers already chosen
    int nb_already_chosen = 1;

    // will hold squared distances to nearest already chosen center
    double *distances = new double[n];

    while (nb_already_chosen < k) {
        double sum_distances = 0.0;
        // calculate squared distance to nearest already chosen center
        for (i = 0; i < n; i++) {
	        distances[i] = DBL_MAX;
            
            for (int j = 0; j < nb_already_chosen; j++) {
                if (points[i].squared_dist(centers[j]) < distances[i])
                    distances[i] = points[i].squared_dist(centers[j]);
            }
            sum_distances += distances[i];
        }

        // choose random point proportional to square distance
        double random = ((double)rand() / RAND_MAX) * sum_distances;
        double prob_sum = 0.0;
        i = 0;
        while (prob_sum <= random && i < n) {
            prob_sum += distances[i];
            i++;
        }
    
        // i-1 is now the index of the chosen point
        for(int m = 0; m < d; m++)
            centers[nb_already_chosen].coords[m] = points[i-1].coords[m];
        nb_already_chosen++;
    }

    delete[] distances;
}

void cloud::init_random_partition() {
    for(int i = 0; i < n; i++) {
        points[i].label = rand() % k;
    }

    set_centroid_centers();
}
