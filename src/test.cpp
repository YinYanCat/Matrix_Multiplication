#include "matrix.cpp"
#include <chrono>
#include <iostream>
int main(){
    std::size_t size = 512;
    Matrix<double> m(size, size);

    for (std::size_t i = 0; i < size; i++) {
        m(i,i) = i;
    }

    auto begin_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::nano> elapsed_time = end_time - begin_time;

    begin_time = std::chrono::high_resolution_clock::now();
    Matrix result =  m.StrassenMultiply(m);
    end_time = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - begin_time);

    std::cout << "Elapsed time: " << elapsed.count() << " ms\n";
    for (std::size_t i = 0; i < size; i++) {
        for (std::size_t j = 0; j < size; j++) {
            if ( i == j ){
                if (result(i,j) != i*i) {
                    std::cout << "Error at " << i << "x" << i << std::endl;
                    break;
                }
            } else {
                if (result(i,j) != 0){
                    std::cout << "Error at " << i << "x" << j << std::endl;
                    break;
                }
            }
        }
    }
    return 0;
}