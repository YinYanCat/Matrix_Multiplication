#pragma once

#include <vector>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <fstream>

bool isPowerofTwo(int n) {    
    return (n > 0) && ((n & (n-1)) == 0);
}

template<typename T>
class Matrix {

    public:
    
    Matrix(std::size_t r, std::size_t c) 
        : cols_(c), rows_(r), data_(c*r, 0) {}
    
    std::size_t cols() const { return cols_; }
    std::size_t rows() const { return rows_; }
    T* data() { return data_.data(); }

    T& operator()(std::size_t r, std::size_t c) {
        return data_[r * cols_ + c];
    }

    T operator()(std::size_t r, std::size_t c) const {
        return data_[r * cols_ + c];
    }

    void saveBin(const std::string& filename) const {
        std::ofstream out(filename, std::ios::binary);
        out.write(reinterpret_cast<const char*>(&rows_), sizeof(rows_));
        out.write(reinterpret_cast<const char*>(&cols_), sizeof(cols_));
        out.write(reinterpret_cast<const char*>(data_.data()), data_.size() * sizeof(T));
    }

    static Matrix loadBin(const std::string& filename) {
        std::ifstream in(filename, std::ios::binary);
        std::size_t r;
        std::size_t c;
        in.read(reinterpret_cast<char*>(&r), sizeof(r));
        in.read(reinterpret_cast<char*>(&c), sizeof(c));
        Matrix m(r,c);
        in.read(reinterpret_cast<char*>(m.data_.data()), m.data_.size() * sizeof(T));

        return m;
    }

    // TODO
    // This should be a const reference, but that would mean that I have to also implement a const view function that makes sures that the view doesn't modify it's matrix.
    Matrix StrassenMultiply(Matrix& a, std::size_t until_size = 0) {
        if (rows_ != cols_ || a.rows() != a.cols()) {
            throw std::runtime_error("One of the matrices isn't square");
        }
        if (cols_ != a.rows()) {
            throw std::runtime_error("The matrices are not the same size");
        }
        if (!isPowerofTwo(rows_)){
            throw std::runtime_error("The matrix size is not a power of 2");
        }
        std::vector<T> Mi(3 * rows_ * cols_, 0);
        Matrix c(rows_, cols_);
        Strassen((*this).view(), a.view(), c.data(), Mi.data(), until_size);
        return c;
    }

    Matrix NaiveMultiply(const Matrix& a){
        Matrix c(rows_, a.cols());
        Naive((*this), a, c.data());
        return c;
    }

    private:
    
    std::size_t cols_;
    std::size_t rows_;
    std::vector<T> data_;

    class MatrixView {

        public:
        MatrixView(T* data, std::size_t r, std::size_t c, std::size_t stride)
            : data_(data), rows_(r), cols_(c), stride_(stride) {}

        std::size_t cols() const { return cols_; }
        std::size_t rows() const { return rows_; }

        T& operator()(std::size_t r, std::size_t c) {
            return data_[r * stride_ + c];
        }

        T operator()(std::size_t r, std::size_t c) const {
            return data_[r * stride_ + c];
        }

        void operator+=(const MatrixView& a) {
            for (std::size_t i = 0; i < rows_; i++) {
                for (std::size_t j = 0; j < cols_; j++) {
                    (*this)(i,j) += a(i,j);
                }
            }
        }

        void operator-=(const MatrixView& a) {
            for (std::size_t i = 0; i < rows_; i++) {
                for (std::size_t j = 0; j < cols_; j++) {
                    (*this)(i,j) -= a(i,j);
                }
            }
        }

        static void flatSum(const MatrixView& a, const MatrixView& b, T* c) {
            std::size_t r_ = a.rows();
            std::size_t c_ = a.cols();
            for (std::size_t i = 0; i < r_; i++) {
                for (std::size_t j = 0; j < c_; j++) {
                    c[i * c_ + j] = a(i,j) + b(i,j);
                }
            }
        }

        static void flatSubs(const MatrixView& a, const MatrixView& b, T* c) {
            std::size_t r_ = a.rows();
            std::size_t c_ = a.cols();
            for (std::size_t i = 0; i < r_; i++) {
                for (std::size_t j = 0; j < c_; j++) {
                    c[i * c_ + j] = a(i,j) - b(i,j);
                }
            }
        }

        void copy(const MatrixView& a) {
            for (std::size_t i = 0; i < rows_; i++) {
                for (std::size_t j = 0; j < cols_; j++) {
                    (*this)(i,j) = a(i,j);
                }
            }
        }

        MatrixView subview(std::size_t rows, std::size_t cols, std::size_t row_offset = 0, std::size_t col_offset = 0) {
            
            return MatrixView(
                    data_ + row_offset * stride_ + col_offset,
                    rows,
                    cols,
                    stride_
            );
        }
        
        static MatrixView flat(T* data, std::size_t rows, std::size_t cols) {
            return MatrixView(data, rows, cols, cols);
        }

        private:
    
        T* data_;
        std::size_t rows_;
        std::size_t cols_;
        std::size_t stride_;

    };

    MatrixView view(std::size_t rows = 0, std::size_t cols = 0, std::size_t row_offset = 0, std::size_t col_offset = 0) {

        if (rows == 0) rows = rows_;
        if (cols == 0) cols = cols_;

        return MatrixView(
                data_.data() + row_offset * cols_ + col_offset,
                rows,
                cols,
                cols_
        );
    }

    template <typename A, typename B>
    static void Naive(const A& a, const B& b, T* c) {
        if (a.cols() != b.rows()) {
            throw std::runtime_error("Invalid operation, the given matrices can't be multiplied.");
        }
        std::size_t a_r = a.rows();
        std::size_t a_c = a.cols();
        std::size_t b_c = b.cols();
        for (std::size_t i = 0; i < a_r; i++) {
            for (std::size_t j = 0; j < b_c; j++) {
                c[i * b_c + j] = 0;
                for (std::size_t k = 0; k < a_c; k++) {
                    c[i * b_c + j] += a(i,k) * b(k,j);
                }
            }
        }
    }

    void Strassen(MatrixView a, MatrixView b, T* c, T* Mi, std::size_t until_size){

        std::size_t n = a.rows();

        // Por motivos de mejor identificación del punto de cruce, preferí cambiar la definición de lo pedido 
        // (Divide & Conquer sólo si n ≥ n0, o Clasico solo si n < n0)
        // a Divide & Conquer solo si n > n0, o Clasico solo si n <= n0
        if (n <= until_size) {  
            Naive(a, b, c);
            return;
        }

        if (a.cols() == 1 && a.rows() == 1 && b.cols() == 1 && b.rows() == 1) {
            c[0] = a(0,0) * b(0,0);
            return;
        }

        std::size_t s = n/2;
        std::size_t s2 = s*s;

        // Punteros a subdivisiones del espacio Mi para calculos
        T* p_a = Mi;
        T* p_b = Mi + s2;
        T* p_M1 = Mi + 2 * s2;
        T* p_M2 = Mi + 3 * s2;
        T* p_M3 = Mi + 4 * s2;
        T* p_M4 = Mi + 5 * s2;
        T* p_M5 = Mi + 6 * s2;
        T* p_M6 = Mi + 7 * s2;
        T* p_M7 = Mi + 8 * s2;
        T* next = Mi + 9 * s2;

        MatrixView A11 = a.subview(s, s);
        MatrixView A12 = a.subview(s, s, 0, s);
        MatrixView A21 = a.subview(s, s, s);
        MatrixView A22 = a.subview(s, s, s, s);

        MatrixView B11 = b.subview(s, s);
        MatrixView B12 = b.subview(s, s, 0, s);
        MatrixView B21 = b.subview(s, s, s);
        MatrixView B22 = b.subview(s, s, s, s);
        
        MatrixView tempA = MatrixView::flat(p_a, s, s);
        MatrixView tempB = MatrixView::flat(p_b, s, s);
        
        MatrixView::flatSum(A11, A22, p_a);
        MatrixView::flatSum(B11, B22, p_b);
        Strassen(tempA, tempB, p_M1, next, until_size);             // (A11 + A22) * (B11 + B22)

        MatrixView::flatSum(A21, A22, p_a);
        Strassen(tempA, B11, p_M2, next, until_size);               // (A21 + A22) * B11

        MatrixView::flatSubs(B12, B22, p_b);
        Strassen(A11, tempB, p_M3, next, until_size);               // A11 * (B12 - B22)

        MatrixView::flatSubs(B21, B11, p_b);
        Strassen(A22, tempB, p_M4, next, until_size);               // A22 * (B21 - B11)

        MatrixView::flatSum(A11, A12, p_a);
        Strassen(tempA, B22, p_M5, next, until_size);               // (A11 + A12) * B22

        MatrixView::flatSubs(A21, A11, p_a);
        MatrixView::flatSum(B11, B12, p_b);
        Strassen(tempA, tempB, p_M6, next, until_size);             // (A21 - A11) * (B11 + B12)

        MatrixView::flatSubs(A12, A22, p_a);
        MatrixView::flatSum(B21, B22, p_b);
        Strassen(tempA, tempB, p_M7, next, until_size);             // (A12 - A22) * (B21 + B22)

        for (std::size_t i = 0; i < s; i++) {
            for (std::size_t j = 0; j < s; j++) {
                std::size_t idx = i * s + j;
                c[i * n + j] = p_M1[idx] + p_M4[idx] - p_M5[idx] + p_M7[idx];               // C11
                c[i * n + (j+s)] = p_M3[idx] + p_M5[idx];                                   // C12
                c[(i+s) * n + j] = p_M2[idx] + p_M4[idx];                                   // C21
                c[(i+s) * n + (j+s)] = p_M1[idx] - p_M2[idx] + p_M3[idx] + p_M6[idx];       // C22
            }
        }        
    }

};