#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "MatrixIndex.hpp"
#include "RBTree.hpp"

class SparseMatrixRowProxy;
class SparseMatrixRow;

class SparseMatrix
{
    static double               eps;
    static size_t               instance_quantity;

    RBTree<MatrixIndex, double> tree;
    size_t                      height; 
    size_t                      width;

    friend SparseMatrixProxy;
    friend SparseMatrixProxy& operator+(size_t, SparseMatrix&);
    friend const SparseMatrixProxy& operator+(size_t, const SparseMatrix&);

    public:
    static void set_precision(double);
    static double get_precision();

    SparseMatrix(size_t, size_t);
    SparseMatrix(const SparseMatrix&);
    ~SparseMatrix();
    SparseMatrix& operator=(const SparseMatrix&);

    double get(size_t, size_t) const;
    void set(size_t, size_t, double);

    size_t num_rows() const;
    size_t num_columns() const;

    SparseMatrixProxy operator[](size_t);
    SparseMatrixProxy operator[](size_t) const;
    SparseMatrixProxy operator*();
    SparseMatrixProxy operator*() const;
    SparseMatrixProxy operator+(size_t);
    SparseMatrixProxy operator+(size_t) const;
    SparseMatrixProxy operator-(size_t);
    SparseMatrixProxy operator-(size_t) const;

    bool operator==(const SparseMatrix&) const;
    bool operator!=(const SparseMatrix&) const;

    /*
    ** Арифметические операции вида SparseMatrix.op(SparseMatrix)
    */
    SparseMatrix& perform_operation(void (*)(double&,double), const SparseMatrix&);
    SparseMatrix& operator+=(const SparseMatrix&);
    SparseMatrix& operator-=(const SparseMatrix&);
    SparseMatrix& operator*=(const SparseMatrix&);
    SparseMatrix& operator/=(const SparseMatrix&);
    SparseMatrix operator+(const SparseMatrix&) const;
    SparseMatrix operator-(const SparseMatrix&) const;
    SparseMatrix operator*(const SparseMatrix&) const;
    SparseMatrix operator/(const SparseMatrix&) const;
    SparseMatrix dot(const SparseMatrix&) const;

    /*
    ** Арифметические операции вида SparseMatrix.op(double)
    */
    SparseMatrix& operator*=(double val);
    SparseMatrix& operator/=(double val);
    SparseMatrix operator*(double val) const;
    SparseMatrix operator/(double val) const;
};

class SparseMatrixProxy {
    bool                        modifyable; // изменяемая ли матрица m
    size_t                      row;
    size_t                      column; 
    unsigned short              deref; // сколько произошло разыменований

    SparseMatrix                *m;
    Node<MatrixIndex,double>    *node;
    size_t                      added; // это страшный костыль

    friend SparseMatrixProxy& operator+(size_t, SparseMatrixProxy&);
    friend c
    public:
    SparseMatrixProxy(SparseMatrix *, size_t, bool = false, bool = true);

    SparseMatrixProxy& operator+(size_t);
    SparseMatrixProxy& operator-(size_t);
    SparseMatrixProxy operator[](size_t);
    SparseMatrixProxy operator*();
    
    SparseMatrixProxy& operator=(double);
    SparseMatrixProxy& operator=(const SparseMatrixRowProxy&);
    
    SparseMatrixProxy& perform_operation(void (*)(double&,double), double);
    SparseMatrixProxy& operator+=(double);
    SparseMatrixProxy& operator-=(double);
    SparseMatrixProxy& operator*=(double);
    SparseMatrixProxy& operator/=(double);

    bool operator==(double) const;
    bool operator==(SparseMatrixProxy&) const;
    bool operator!=(double) const;
    bool operator!=(SparseMatrixProxy&) const;

    operator double() const;
}

#endif