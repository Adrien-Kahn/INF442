#include "point.hpp"
#include "cloud.hpp"
#include "kernel.hpp"
#include "knn.hpp"
#include <iostream>

#include <cfloat>


knn::knn(cloud *data, int k_, double V_) : kernel(data) {

k = k_;
V = V_;

}


double knn::density(point &p) {
	
	int minidx = -1;
	double min = DBL_MAX;
	int n = data->get_n();
	point *points = data->points;
	
	int idx[k];
	for (int i = 0; i < k; i++) idx[i] = -1;
	
	for (int i = 0; i < k; i++) {
		
		int minidx = -1;
		double min = DBL_MAX;
		
		for (int j = 0; j < n; j++) {
			
			bool b = true;
			
			for (int ii = 0; ii < k; ii++) {
				if (j == idx[ii]) b = false;
			}
			
			if (points[j].dist(p) < min && b) {
				minidx = j;
				min = points[j].dist(p);
			}
		
		}
		
		idx[i] = minidx;
		
	}
	
	return k / (2 * n * V * points[idx[k - 1]].dist(p));

}
