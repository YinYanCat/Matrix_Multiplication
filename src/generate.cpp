#include "../src/matrix.hpp"
#include <random>
#include <string>
#include <fstream>

Matrix<int> genIntMatrix(std::size_t n){
    Matrix<int> result(n,n);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(-1000, 1000);
    for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < n; j++){
            result(i,j) = dist(gen);
        }
    }
    return result;
}

Matrix<double> genDoubleMatrix(std::size_t n){
    Matrix<double> result(n,n);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < n; j++){
            result(i,j) = dist(gen);
        }
    }
    return result;
}

Matrix<int> genDiagonalIntMatrix(std::size_t n){
    Matrix<int> result(n,n);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(-1000, 1000);
    for (std::size_t i = 0; i < n; i++) {
        result(i,i) = dist(gen);
    }
    return result;
}

Matrix<double> genDiagonalDoubleMatrix(std::size_t n){
    Matrix<double> result(n,n);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (std::size_t i = 0; i < n; i++) {
        result(i,i) = dist(gen);
    }
    return result;
}

void genSaveTests(std::size_t n, int quantity) {
    for (int i = 0; i < quantity; i++) {
        std::string path = std::to_string(n) + "/" + std::to_string(i);

        Matrix<int> r_i  = genIntMatrix(n);
        Matrix<double> r_d  = genDoubleMatrix(n);
        Matrix<int> r_di = genDiagonalIntMatrix(n);
        Matrix<double> r_dd = genDiagonalDoubleMatrix(n);

        r_i .saveBin("int/" + path);
        r_d .saveBin("double/" + path);
        r_di.saveBin("diagonal_int/" + path);
        r_dd.saveBin("diagonal_double/" + path);

    }
}
