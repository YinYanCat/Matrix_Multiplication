#pragma once

#include "./matrix.cpp"


Matrix Strassen(MatrixView a, MatrixView b){
    
    if (a.cols() != b.rows()) {
        return Matrix(1,1);
    }

    if (a.cols() == 1 && a.rows() == 1 && b.cols() == 1 && b.rows() == 1) {
        Matrix c(1,1);
        c(0,0) = a(0,0) * b(0,0);
        return c;
    }

    Matrix c = Matrix(a.rows(), b.cols());
    
    Matrix M1 = Strassen(
        (
            a.subview(a.rows()/2, a.cols()/2) +                             //A11
            a.subview(a.rows()/2, a.cols()/2, a.rows()/2, a.cols()/2)       //A22
        ).view(),
        (
            b.subview(b.rows()/2, b.cols()/2) +                             //B11
            b.subview(b.rows()/2, b.cols()/2, b.rows()/2, b.cols()/2)       //B22
        ).view()
    );

    Matrix M2 = Strassen(
        (
            a.subview(a.rows()/2, a.cols()/2, a.rows()/2) +                 //A21
            a.subview(a.rows()/2, a.cols()/2, a.rows()/2, a.cols()/2)       //A22
        ).view(), 
        (
            b.subview(b.rows()/2, b.cols()/2)
        )
    );

    Matrix M3 = Strassen(
        (
            a.subview(a.rows()/2, a.cols()/2)                               //A11
            
        ),
        (
            b.subview(b.rows()/2, b.cols()/2, 0, b.cols()/2) -              //B12
            b.subview(b.rows()/2, b.cols()/2, b.rows()/2, b.cols()/2)       //B22
        ).view()
    );

    Matrix M4 = Strassen(
        (
            a.subview(a.rows()/2, a.cols()/2, a.rows()/2, a.cols()/2)       //A22
        ),
        (
            b.subview(b.rows()/2, b.cols()/2, b.rows()/2) -                 //B21
            b.subview(b.rows()/2, b.cols()/2)                               //B11
        ).view()
    );

    Matrix M5 = Strassen(
        (
            a.subview(a.rows()/2, a.cols()/2) +                             //A11
            a.subview(a.rows()/2, a.cols()/2, 0, a.cols()/2)                //A12
        ).view(),
        (
            b.subview(b.rows()/2, b.cols()/2, b.rows()/2, b.cols()/2)       //B22
        )
    );
    Matrix M6 = Strassen(
        (
            a.subview(a.rows()/2, a.cols()/2, a.rows()/2) -                 //A21
            a.subview(a.rows()/2, a.cols()/2)                               //A11
        ).view(),
        (
            b.subview(b.rows()/2, b.cols()/2) +                             //B11
            b.subview(b.rows()/2, b.cols()/2, 0, b.cols()/2)                //B12
        ).view()
    );
    Matrix M7 = Strassen(
        (
            a.subview(a.rows()/2, a.cols()/2, 0, a.cols()/2) -              //A12
            a.subview(a.rows()/2, a.cols()/2, a.rows()/2, a.cols()/2)       //A22
        ).view(),
        (
            b.subview(b.rows()/2, b.cols()/2, b.rows()/2) +                 //B21
            b.subview(b.rows()/2, b.cols()/2, b.rows()/2, b.cols()/2)       //B22
        ).view()
    );

    MatrixView C11 = c.view(c.rows()/2, c.cols()/2);
    MatrixView C12 = c.view(c.rows()/2, c.cols()/2, 0, c.cols()/2);
    MatrixView C21 = c.view(c.rows()/2, c.cols()/2, c.rows()/2);
    MatrixView C22 = c.view(c.rows()/2, c.cols()/2, c.rows()/2, c.cols()/2);

    C11.assign(M1 + M4 + M7 - M5);
    C12.assign(M3 + M5);
    C21.assign(M2 + M4);
    C22.assign(M1 - M2 + M3 + M6);
    
    return c;
}

Matrix StrassenMultiply(Matrix& a, Matrix& b){
    return Strassen(a.view(), b.view());
}