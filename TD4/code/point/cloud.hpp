#pragma once  // ensures that this header file is only included once

#include "point.hpp"  // include the header file with the declaration
		      // of the class point, so that it can be used
		      // here

class cloud	      // this is just a declaration of the class
		      // implementation is in a separate file:
		      // cloud.cpp
{
private:
  int d;
  int n;
  int k;

  // maximum possible number of points
  int nmax;

  point *points;
  point *centers;

public:
  cloud(int _d, int _nmax, int _k);
  ~cloud();

  // Getters
  
  int get_d();
  int get_n();
  int get_k();

  point& get_point(int i);
  point& get_center(int j);

  // Setters
  
  void set_center(point &p, int j);

  // Helper methods
  
  void add_point(point &p, int label);

  double intracluster_variance();

  int set_voronoi_labels();
  void set_centroid_centers();

  void kmeans();

  // Cluster initialisation options
  
  void init_forgy();
  void init_plusplus();
  void init_random_partition();
};
