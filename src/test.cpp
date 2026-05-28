#include "matrix.hpp"
#include <chrono>
#include <iostream>
int main(){
    std::size_t size = 4;
    Matrix<int> m(size, size);


    for (std::size_t i = 0; i < size; i++) {
        for (std::size_t j = 0; j < size; j++) {
            m(i,j) = 1;
        }
        std::cout << std::endl;
    }

    auto begin_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::nano> elapsed_time = end_time - begin_time;

    begin_time = std::chrono::high_resolution_clock::now();
    Matrix result =  m.NaiveMultiply(m); //Strassen puro
    end_time = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time);

    std::cout << "input matrix M:" << std::endl;
    for (std::size_t i = 0; i < size; i++) {
        for (std::size_t j = 0; j < size; j++) {
            std::cout << m(i,j) << "\t";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "M x M:" << std::endl;
    for (std::size_t i = 0; i < size; i++) {
        for (std::size_t j = 0; j < size; j++) {
            std::cout << result(i,j) << "\t";
        }
        std::cout << std::endl;
    }
    return 0;
}