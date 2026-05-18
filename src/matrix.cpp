#pragma once

#include <vector>
#include <cstddef>


class MatrixView;

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

    Matrix& operator+=(const Matrix& a)  {
        for (std::size_t i = 0; i<data_.size(); i++){
            data_[i] += a.data_[i];
        }
        return *this;
    }

    Matrix& operator-=(const Matrix& a)  {
        for (std::size_t i = 0; i<data_.size(); i++){
            data_[i] -= a.data_[i];
        }
        return *this;
    }

    Matrix operator+(const Matrix& a) const {
        Matrix result = *this;
        result += a;
        return result;
    }

    Matrix operator-(const Matrix& a) const {
        Matrix result = *this;
        result -= a;
        return result;
    }

    MatrixView view(std::size_t rows = 0, std::size_t cols = 0, std::size_t row_offset = 0, std::size_t col_offset = 0);

    private:
    
    std::size_t cols_;
    std::size_t rows_;
    std::vector<double> data_;

};


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

    Matrix operator+(const MatrixView& a) const{
        
        Matrix result(rows_, cols_);

        for(std::size_t i = 0; i < rows_; i++){
            for(std::size_t j = 0; j < cols_; j++){
                result(i,j) = (*this)(i,j) + a(i,j);
            }
        }

        return result;
    }

    Matrix operator-(const MatrixView& a) const{
        
        Matrix result(rows_, cols_);

        for(std::size_t i = 0; i < rows_; i++){
            for(std::size_t j = 0; j < cols_; j++){
                result(i,j) = (*this)(i,j) - a(i,j);
            }
        }

        return result;
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

    MatrixView Matrix::view(std::size_t rows, std::size_t cols, std::size_t row_offset, std::size_t col_offset) {

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