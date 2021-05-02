#include "Classification.hpp"
#include "Dataset.hpp"

Classification::Classification(Dataset* dataset, int col_class) {
    m_dataset = dataset;
    m_col_class = col_class;
}

Dataset* Classification::getDataset(){
    return m_dataset;
}

int Classification::getColClass(){
    return m_col_class;
}
