#pragma once

#include "./matrix.cpp"

Matrix NaiveMultiply(Matrix a, Matrix b){
    
    if (a.cols() != b.rows()) {
        return;
    }

    Matrix result = Matrix(a.rows(), b.cols());

    for (std::size_t i = 0; i < a.rows(); i++) {
        for (std::size_t j = 0; j < b.cols(); j++) {
            for (std::size_t k = 0; k < a.cols(); k++) {
                result(i,j) += a(i,k) * b(k,j);
            }
        }
    }
    
    return result;
}