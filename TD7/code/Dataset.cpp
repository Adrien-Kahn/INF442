
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include "Dataset.hpp"


int Dataset::GetNbrSamples() const {
	return m_nsamples;
}

int Dataset::GetDim() const {
	return m_dim;
}

Dataset::~Dataset() {
}

void Dataset::Show(bool verbose) const {
	std::cout<<"Dataset with "<<m_nsamples<<" samples, and "<<m_dim<<" dimensions."<<std::endl;

	if (verbose) {
		for (int i=0; i<m_nsamples; i++) {
			for (int j=0; j<m_dim; j++) {
				std::cout<<m_instances[i][j]<<" ";
			}
			std::cout<<std::endl;		
		}
	}
}

Dataset::Dataset(const char* file) {
	m_nsamples = 0;
	m_dim = -1;

	std::ifstream fin(file);
	
	if (fin.fail()) {
		std::cout<<"Cannot read from file "<<file<<" !"<<std::endl;
		exit(1);
	}
	
	std::vector<double> row; 
    std::string line, word, temp; 

	while (getline(fin, line)) {
		row.clear();
        std::stringstream s(line);
        
        int valid_sample = 1;
        int ncols = 0;
        while (getline(s, word, ',')) { 
            // add all the column data 
            // of a row to a vector 
            double val = std::atof(word.c_str());
            row.push_back(val);
            ncols++;
        }
        if (!valid_sample) {
            // in this version, valid_sample is always 1
        	continue;
        }         
        m_instances.push_back(row);
        if (m_dim==-1) m_dim = ncols;
        else if (m_dim!=ncols) {
        	std::cerr << "ERROR, inconsistent dataset" << std::endl;
        	exit(-1);
        }
        
		m_nsamples ++;
	}
	
	fin.close();
}

const std::vector<double>& Dataset::GetInstance(int i) const {
	return m_instances[i];
}
