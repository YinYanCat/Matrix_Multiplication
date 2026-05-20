#include <vector>
#include <cstddef>
#include <stdexcept>

class Matrix;

template <typename A, typename B, typename Oprt>
Matrix add(const A& a, const B& b, Oprt operation);

class Matrix {

    public:
    Matrix(std::size_t r, std::size_t c) 
        : cols_(c), rows_(r), data_(c*r) {}
    
    std::size_t cols() const { return cols_; }
    std::size_t rows() const { return rows_; }

    double& operator()(std::size_t r, std::size_t c) {
        return data_[r * cols_ + c];
    }

    double operator()(std::size_t r, std::size_t c) const {
        return data_[r * cols_ + c];
    }

    Matrix operator+(const Matrix& a) const {
        return add((*this), a, std::plus<double>());
    }

    Matrix operator-(const Matrix& a) const {
        return add((*this), a, std::minus<double>());
    }

    // TODO
    // This should be a const reference, but that would mean that I have to also implement a const view function that makes sures that the view doesn't modify it's matrix.
    Matrix StrassenMultiply(Matrix& a){
        return Strassen((*this).view(), a.view());
    }

    Matrix NaiveMultiply(const Matrix& a){
    
        if (cols_ != a.rows()) {
            throw std::runtime_error("Invalid operation, the given matrices can't be multiplied.");
        }

        Matrix result = Matrix(rows_, a.cols());

        for (std::size_t i = 0; i < rows_; i++) {
            for (std::size_t j = 0; j < a.cols(); j++) {
                for (std::size_t k = 0; k < cols_; k++) {
                    result(i,j) += (*this)(i,k) * a(k,j);
                }
            }
        }
        
        return result;
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

        Matrix operator+(const MatrixView& a) const {
            return add(*this, a, std::plus<double>());
        }
        Matrix operator-(const MatrixView& a) const {
            return add(*this, a, std::minus<double>());
        }

        void assign(const Matrix& a) {
            for(std::size_t i = 0; i < rows_; i++) {
                for(std::size_t j = 0; j < cols_; j++) {
                    (*this)(i,j) = a(i,j);
                }
            }
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

    Matrix Strassen(MatrixView a, MatrixView b){
    
        if (a.cols() != b.rows()) {
            throw std::runtime_error("Invalid operation, the given matrices can't be multiplied.");
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

};

template <typename A, typename B, typename Oprt> 
Matrix add(const A &a,const B &b, Oprt operation){
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