#include "SparseMatrix.hpp"
#include <iostream>
#include <stdexcept>

int main() {
    SparseMatrix::set_precision(1e-6);
    SparseMatrix m(5, 5);
    
    for (size_t i = 0; i < m.num_rows(); i++) {
        for (size_t j = 0; j < m.num_columns(); j++) {
            *((m + i)[0] + j) = ((*(m + j))[i] = i + j);
        }
    }

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            std::cout << m[i][j] << ' ';
        }
        std::cout << '\n';
    }
    std::cout << '\n';
    
    const SparseMatrix mm(m);
    try {
        for (size_t i = 0; i < mm.num_columns(); i++) {
            for (size_t j = 0; j < mm.num_rows(); j++) {
                *(mm[i] + j) += (*(mm[i] + j) + 2 - 2);
                *(mm[i] + j) /= 2;
                *(mm[i] + j) *= (2 + *(mm[i] + j) - 2);
            }
            std::cout << '\n';
        }
    } catch (const std::runtime_error& e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    try {
        SparseMatrix::set_precision(1e-10);
    } catch (const std::runtime_error& e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    try {
        m[5][4];
    } catch (const std::out_of_range& e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
    try {
        m[4][5];
    } catch (const std::out_of_range& e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }

    SparseMatrix mmm(1,1);
    mmm = mm;
    for (size_t i = 0; i < mmm.num_rows(); i++) {
        std::cout << ((1 + mmm[i] + i - 2)[1] -= 9) << ' ';
    }
    std::cout << "\n\n";

    mmm = (mmm *= 2) * (mmm  /= 2);
    for (size_t i = 0; i < mmm.num_columns(); i++) {
        for (size_t j = 0; j < mmm.num_rows(); j++) {
            std::cout << (*(mmm + i))[j] << ' ';
        }
        std::cout << '\n';
    }
    std::cout << '\n';
    m = m.dot(m);
    for (size_t i = 0; i < m.num_columns(); i++) {
        for (size_t j = 0; j < m.num_rows(); j++) {
            std::cout << (*(m + i))[j] << ' ';
        }
        std::cout << '\n';
    }
    std::cout << '\n';
    std::cout << (m * 2 == m + m ? "true" : "false") << std::endl;
    m[1][1] = m[2][2] = m[3][3] = m[4][4] = 0.0;
    m[1][1] = m[2][2] = m[3][3] = m[4][4] = 10000.0;
    for (size_t i = 0; i < m.num_columns(); i++) {
        for (size_t j = 0; j < m.num_rows(); j++) {
            std::cout << (*(m + i))[j] << ' ';
        }
        std::cout << '\n';
    }
    m = m;
    std::cout << '\n';
    std::cout << **m << '\n';
    // m.set(5, 1, 1.0);
    return 0;
}