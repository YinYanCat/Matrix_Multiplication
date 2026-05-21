#include <vector>
#include <cstddef>
#include <stdexcept>

class Matrix {

    public:
    
    Matrix(std::size_t r, std::size_t c) 
        : cols_(c), rows_(r), data_(c*r,0) {}
    
    std::size_t cols() const { return cols_; }
    std::size_t rows() const { return rows_; }

    double& operator()(std::size_t r, std::size_t c) {
        return data_[r * cols_ + c];
    }

    double operator()(std::size_t r, std::size_t c) const {
        return data_[r * cols_ + c];
    }

    // TODO
    // This should be a const reference, but that would mean that I have to also implement a const view function that makes sures that the view doesn't modify it's matrix.
    Matrix StrassenMultiply(Matrix& a, std::size_t until_size = 0){
        return Strassen((*this).view(), a.view(), until_size);
    }

    Matrix NaiveMultiply(const Matrix& a){
        return Naive((*this), a);
    }

    private:
    
    std::size_t cols_;
    std::size_t rows_;
    std::vector<double> data_;

    class MatrixView {

        public:
        MatrixView(Matrix* m, std::size_t r, std::size_t c, std::size_t ro = 0, std::size_t co = 0)
            : m(m), rows_(r), cols_(c), row_offset(ro), col_offset(co) {}

        std::size_t cols() const { return cols_; }
        std::size_t rows() const { return rows_; }

        double& operator()(std::size_t r, std::size_t c) {
            return (*m)(row_offset + r, col_offset + c);
        }

        double operator()(std::size_t r, std::size_t c) const {
            return (*m)(row_offset + r, col_offset + c);
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

        Matrix operator+(const MatrixView& a) const { return add(*this, a, std::plus<double>()); }
        Matrix operator-(const MatrixView& a) const { return add(*this, a, std::minus<double>()); }

        void copy(const MatrixView& a) {
            for (std::size_t i = 0; i < rows_; i++) {
                for (std::size_t j = 0; j < cols_; j++) {
                    (*this)(i,j) = a(i,j);
                }
            }
        }

        Matrix NaiveMultiply(const MatrixView& a){
            return Naive((*this), a);
        }

        MatrixView subview(std::size_t rows = 0, std::size_t cols = 0, std::size_t row_offset = 0, std::size_t col_offset = 0) {

            if (rows == 0) rows = rows_;
            if (cols == 0) cols = cols_;

            return MatrixView(
                    m,
                    rows,
                    cols,
                    this->row_offset + row_offset,
                    this->col_offset + col_offset
            );
        }
        
        private:
        
        Matrix* m;
        std::size_t rows_;
        std::size_t cols_;
        std::size_t row_offset;
        std::size_t col_offset;

        template <typename Oprt> 
        static Matrix add(const MatrixView& a, const MatrixView& b, Oprt operation){
            std::size_t rows_ = a.rows();
            std::size_t cols_ = a.cols();

            if (rows_ != b.rows() || cols_ != b.cols()) {
                throw std::runtime_error("Invalid operation, the given matrices have to be the same size.");
            }
            Matrix result(rows_, cols_);
            for(std::size_t i = 0; i < rows_; i++){
                for(std::size_t j = 0; j < cols_; j++){
                    result(i,j) = operation(a(i,j),b(i,j));
                }
            }
            return result;
        }

    };

    MatrixView view(std::size_t rows = 0, std::size_t cols = 0, std::size_t row_offset = 0, std::size_t col_offset = 0) {

        if (rows == 0) rows = rows_;
        if (cols == 0) cols = cols_;

        return MatrixView(
                this,
                rows,
                cols,
                row_offset,
                col_offset
        );
    }

    template <typename A, typename B>
    static Matrix Naive(const A& a, const B& b) {
        if (a.cols() != b.rows()) {
            throw std::runtime_error("Invalid operation, the given matrices can't be multiplied.");
        }
        Matrix c(a.rows(), b.cols());
        for (std::size_t i = 0; i < a.rows(); i++) {
            for (std::size_t j = 0; j < b.cols(); j++) {
                for (std::size_t k = 0; k < a.cols(); k++) {
                    c(i,j) += a(i,k) * b(k,j);
                }
            }
        }
        return c;
    }

    Matrix Strassen(MatrixView a, MatrixView b, std::size_t until_size){
    
        if (a.cols() != b.rows()) {
            throw std::runtime_error("Invalid operation, the given matrices can't be multiplied.");
        }

        // TODO
        // Implement either matrix padding or checking for square matrices of size power 2

        // Here it checks if the row size is equal to the stop size to start using the naive algorithm
        // At the moment of implementing this there is no checking on the sizes of the matrices to be
        // a power of 2

        if (a.rows() < until_size) {
            return a.NaiveMultiply(b);
        }

        if (a.cols() == 1 && a.rows() == 1 && b.cols() == 1 && b.rows() == 1) {
            Matrix c(1,1);
            c(0,0) = a(0,0) * b(0,0);
            return c;
        }

        Matrix c = Matrix(a.rows(), b.cols());

        std::size_t ahr = a.rows()/2;                                                       // a half rows
        std::size_t ahc = a.cols()/2;                                                       // a half columns
        std::size_t bhr = b.rows()/2;                                                       // b half rows
        std::size_t bhc = b.cols()/2;                                                       // b half columns

        MatrixView A11 = a.subview(ahr, ahc);
        MatrixView A12 = a.subview(ahr, ahc, 0, ahc);
        MatrixView A21 = a.subview(ahr, ahc, ahr, 0);
        MatrixView A22 = a.subview(ahr, ahc, ahr, ahc);

        MatrixView B11 = b.subview(bhr, bhc);
        MatrixView B12 = b.subview(bhr, bhc, 0, bhc);
        MatrixView B21 = b.subview(bhr, bhc, bhr, 0);
        MatrixView B22 = b.subview(bhr, bhc, bhr, bhc);

        MatrixView C11 = c.view(c.rows()/2, c.cols()/2);
        MatrixView C12 = c.view(c.rows()/2, c.cols()/2, 0, c.cols()/2);
        MatrixView C21 = c.view(c.rows()/2, c.cols()/2, c.rows()/2);
        MatrixView C22 = c.view(c.rows()/2, c.cols()/2, c.rows()/2, c.cols()/2);
        
        Matrix temp1 = A11 + A22;
        Matrix temp2 = B11 + B22;
        Matrix M1 = Strassen(temp1.view(), temp2.view(), until_size);       // (A11 + A22) * (B11 + B22)

        temp1 = A21 + A22;
        Matrix M2 = Strassen(temp1.view(), B11, until_size);                // (A21 + A22) * B11

        temp2 = B12 - B22;
        Matrix M3 = Strassen(A11, temp2.view(), until_size);                // A11 * (B12 - B22)

        temp2 = B21 - B11;
        Matrix M4 = Strassen(A22, temp2.view(), until_size);                // A22 * (B21 - B11)

        temp1 = A11 + A12;
        Matrix M5 = Strassen(temp1.view(), B22, until_size);                // (A11 + A12) * B22

        temp1 = A21 - A11;
        temp2 = B11 + B12;
        Matrix M6 = Strassen(temp1.view(), temp2.view(), until_size);       // (A21 - A11) * (B11 + B12)

        temp1 = A12 - A22;
        temp2 = B21 + B22;
        Matrix M7 = Strassen(temp1.view(), temp2.view(), until_size);       // (A12 - A22) * (B21 + B22)

        C11.copy(M1.view());
        C11 += M4.view();
        C11 += M7.view();
        C11 -= M5.view();

        C12.copy(M3.view());
        C12 += M5.view();

        C21.copy(M2.view());
        C21 += M4.view();

        C22.copy(M1.view());
        C22 -= M2.view();
        C22 += M3.view();
        C22 += M6.view();
        
        return c;
    }

};