#include "SparseMatrix.hpp"
#include <iostream>

int main() {
    SparseMatrix m(5, 5);

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            std::cout << m[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    for (size_t i = 0; i < m.num_rows(); i++) {
        for (size_t j = 0; j < m.num_columns(); j++) {
            *((m + i - 2)[2] + j) = ((*(m + j))[i] = i + j);
        }
    }
    
    const SparseMatrix mm(m);
    for (size_t i = 0; i < mm.num_columns(); i++) {
        for (size_t j = 0; j < mm.num_rows(); j++) {
            std::cout << *(mm[i] + j) + 10 << ' ';
            *(mm[i] + j) += *(mm[i] + j);
            *(mm[i] + j) /= 2;
            *(mm[i] + j) *= *(mm[i] + j);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << SparseMatrix::get_precision() << std::endl;
    SparseMatrix::set_precision(1e-10);
    std::cout << SparseMatrix::get_precision() << std::endl;
    std::cout << mm.num_columns() << std::endl;
    std::cout << mm.num_rows() << std::endl << std::endl;

    SparseMatrix mmm(1,1);
    mmm = mm;
    for (size_t i = 0; i < mmm.num_rows(); i++) {
        (mmm[i] + i - 1)[1] -= 9;
    }

    for (size_t i = 0; i < mmm.num_columns(); i++) {
        for (size_t j = 0; j < mmm.num_rows(); j++) {
            std::cout << (*(mmm + i))[j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    mmm *= 2;
    for (size_t i = 0; i < mmm.num_columns(); i++) {
        for (size_t j = 0; j < mmm.num_rows(); j++) {
            std::cout << (*(mmm + i))[j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    mmm = (mmm / 2) * (mmm  / 2);
    for (size_t i = 0; i < mmm.num_columns(); i++) {
        for (size_t j = 0; j < mmm.num_rows(); j++) {
            std::cout << (*(mmm + i))[j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    m = m.dot(m);
    for (size_t i = 0; i < m.num_columns(); i++) {
        for (size_t j = 0; j < m.num_rows(); j++) {
            std::cout << (*(m + i))[j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << (m * 2 == m + m) << std::endl;
    m[1][1] = m[2][2] = m[3][3] = m[4][4] = 0.0;
    m[1][1] = m[2][2] = m[3][3] = m[4][4] = 10000.0;
    for (size_t i = 0; i < m.num_columns(); i++) {
        for (size_t j = 0; j < m.num_rows(); j++) {
            std::cout << (*(m + i))[j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    return 0;
}