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
    size_t                      row;

    friend SparseMatrixRow;
    friend SparseMatrixRowProxy;

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

    SparseMatrixRow operator[](size_t);
    SparseMatrixRow operator[](size_t) const;
    SparseMatrixRow operator*();
    SparseMatrixRow operator*() const;
    SparseMatrix& operator+(size_t);
    const SparseMatrix& operator+(size_t) const;
    SparseMatrix& operator-(size_t);
    const SparseMatrix& operator-(size_t) const;

    bool operator==(const SparseMatrix&) const; // after indexation
    bool operator!=(const SparseMatrix&) const; // after indexation

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

class SparseMatrixRow
{
    bool                        modifyable;
    size_t                      row;
    size_t                      column;
    SparseMatrix                *m;
    Node<MatrixIndex,double>    *node;

    friend SparseMatrixRowProxy;
    friend SparseMatrixRow SparseMatrixRow::operator[](size_t);
    friend SparseMatrixRow SparseMatrixRow::operator[](size_t) const;
    friend SparseMatrixRow SparseMatrixRow::operator*();
    friend SparseMatrixRow SparseMatrixRow::operator*() const;

    SparseMatrixRow(SparseMatrix *, size_t, bool = true);
    
    public:
    SparseMatrixRow& operator+(size_t);
    SparseMatrixRow& operator-(size_t);
    SparseMatrixRowProxy operator[](size_t);
    SparseMatrixRowProxy operator*();
};

class SparseMatrixRowProxy
{
    SparseMatrixRow   *real;

    friend SparseMatrixRowProxy SparseMatrixRowProxy::operator[](size_t);
    friend SparseMatrixRowProxy SparseMatrixRowProxy::operator*();

    SparseMatrixRowProxy(SparseMatrixRow *);

    public:
    SparseMatrixRowProxy& operator=(double);
    SparseMatrixRowProxy& operator=(const SparseMatrixRowProxy&);
    
    SparseMatrixRowProxy& perform_operation(void (*)(double&,double), double);
    SparseMatrixRowProxy& operator+=(double);
    SparseMatrixRowProxy& operator-=(double);
    SparseMatrixRowProxy& operator*=(double);
    SparseMatrixRowProxy& operator/=(double);

    bool operator==(double) const;
    bool operator==(SparseMatrixRowProxy&) const;
    bool operator!=(double) const;
    bool operator!=(SparseMatrixRowProxy&) const;

    operator double() const;
};

#endif