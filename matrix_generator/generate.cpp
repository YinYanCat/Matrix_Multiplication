#include "../src/matrix.cpp"
#include <random>
#include <string>
#include <fstream>

void genMatrix(std::size_t n, int quantity) {
    for (int i = 0; i < quantity; i++) {
        Matrix<double> result(n,n);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        // diagonal matrix
        for (std::size_t j = 0; j < n; j++){
            result(j,j) = dist(gen);
        }
        result.saveBin("output_matrix/diagonal_"+std::to_string(i)+"_"+std::to_string(n)+"x"+std::to_string(n));

        Matrix<double> temp = result;

        // upper triangular
        for (std::size_t j = 0; j < n; j++){
            for (std::size_t k = j + 1; k < n; k++){
                result(j,k) = dist(gen);
            }
        }
        result.saveBin("output_matrix/uptriangular_"+std::to_string(i)+"_"+std::to_string(n)+"x"+std::to_string(n));

        //complete
        for (std::size_t j = 0; j < n; j++){
            for (std::size_t k = j + 1; k < n; k++){
                result(k,j) = dist(gen);
            }
        }
        result.saveBin("output_matrix/full_"+std::to_string(i)+"_"+std::to_string(n)+"x"+std::to_string(n));

        result = temp;
        //lower triangular
        for (std::size_t j = 0; j < n; j++){
            for (std::size_t k = j + 1; k < n; k++){
                result(k,j) = dist(gen);
            }
        }
        result.saveBin("output_matrix/lotriangular_"+std::to_string(i)+"_"+std::to_string(n)+"x"+std::to_string(n));
    }
}

int main(){
    genMatrix(2048,3);
    return 0;
}