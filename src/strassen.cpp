#pragma once

#include "./matrix.cpp"

Matrix StrassenMultiply(Matrix a, Matrix b){
    Strassen(a.view(), b.view());
}

Matrix Strassen(MatrixView a, MatrixView b){
    
    if (a.cols() != b.rows()) {
        return;
    }

    if (a.cols() == 1 && a.rows() == 1 && b.cols() == 1 && b.rows() == 1) {
        Matrix result(1,1);
        result(0,0) = a(0,0) * b(0,0);
        return result;
    }

    Matrix result = Matrix(a.rows(), b.cols());

    Matrix M1 = Matrix(a.rows()/2, b.cols()/2);
    Matrix M2 = Matrix(a.rows()/2, b.cols()/2);
    Matrix M3 = Matrix(a.rows()/2, b.cols()/2);
    Matrix M4 = Matrix(a.rows()/2, b.cols()/2);
    Matrix M5 = Matrix(a.rows()/2, b.cols()/2);
    Matrix M6 = Matrix(a.rows()/2, b.cols()/2);
    Matrix M7 = Matrix(a.rows()/2, b.cols()/2);
    
    M1 = StrassenMultiply((a.subview(a.rows()/2, a.cols()/2) + a.subview(a.rows()/2, a.cols()/2, a.rows()/2, a.cols()/2)), b.subview(b.rows()/2, b.cols()/2) + b.subview(b.rows()/2, b.cols()/2, b.rows()/2, b.cols()/2));
    
    for (std::size_t i = 0; i < a.rows(); i++) {
        for (std::size_t j = 0; j < b.cols(); j++) {
            for (std::size_t k = 0; k < a.cols(); k++) {
                result(i,j) += a(i,k) * b(k,j);
            }
        }
    }
    
    return result;
}