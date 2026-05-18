#pragma once

#include "./matrix.cpp"
#include <stdexcept>

Matrix Strassen(MatrixView a, MatrixView b){
    
    if (a.cols() != b.rows()) {
        throw std::runtime_error("Invalid operation, the given matrices can't be multiplied");
    }

    if (a.cols() == 1 && a.rows() == 1 && b.cols() == 1 && b.rows() == 1) {
        Matrix c(1,1);
        c(0,0) = a(0,0) * b(0,0);
        return c;
    }

    Matrix c = Matrix(a.rows(), b.cols());

    std::size_t ahr = a.rows()/2;                                           //a half rows
    std::size_t ahc = a.cols()/2;                                           //a half columns
    std::size_t bhr = b.rows()/2;                                           //b half rows
    std::size_t bhc = b.cols()/2;                                           //b hals columns

    
    Matrix temp1 = a.subview(ahr, ahc) +                                    //A11 +
            a.subview(ahr, ahc, ahr, ahc);                                  //A22

    Matrix temp2 = b.subview(bhr, bhc) +                                    //B11 +
            b.subview(bhr, bhc, bhr, bhc);                                  //B22

    Matrix M1 = Strassen(temp1.view(), temp2.view());                       //(A11 + A22)*(B11 + B22)

    temp1 = a.subview(ahr, ahc, ahr) +                                      //A21 +
            a.subview(ahr, ahc, ahr, ahc);                                  //A22

    Matrix M2 = Strassen(temp1.view(), b.subview(bhr, bhc));                //(A21 + A22) * B11

    temp2 = b.subview(bhr, bhc, 0, bhc) -                                   //B12 -
            b.subview(bhr, bhc, bhr, bhc);                                  //B22

    Matrix M3 = Strassen(a.subview(ahr, ahc), temp2.view());                //A11 * (B12 - B22)

    temp2 = b.subview(bhr, bhc, bhr) -                                      //B21 -
            b.subview(bhr, bhc);                                            //B11

    Matrix M4 = Strassen(a.subview(ahr, ahc, ahr, ahc), temp2.view());      //A22 * (B12 - B22)

    temp1 = a.subview(ahr, ahc) +                                           //A11 +
            a.subview(ahr, ahc, 0, ahc);                                    //A12

    Matrix M5 = Strassen(temp1.view(), b.subview(bhr, bhc, bhr, bhc));      //(A11 + A12) * B22

    temp1 = a.subview(ahr, ahc, ahr) -                                      //A21 -
            a.subview(ahr, ahc);                                            //A11

    temp2 = b.subview(bhr, bhc) +                                           //B11 +
            b.subview(bhr, bhc, 0, bhc);                                    //B12

    Matrix M6 = Strassen(temp1.view(), temp2.view());                       //(A21 - A11) * (B11 + B12)

    temp1 = a.subview(ahr, ahc, 0, ahc) -                                   //A12 -
            a.subview(ahr, ahc, ahr, ahc);                                  //A22

    temp2 = b.subview(bhr, bhc, bhr) +                                      //B21 +
            b.subview(bhr, bhc, bhr, bhc);                                  //B22

    Matrix M7 = Strassen(temp1.view(), temp2.view());                       //(A12 - A22) * (B21 + B22)

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