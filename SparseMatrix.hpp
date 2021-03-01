#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "MatrixIndex.hpp"
#include "RBTree.hpp"

class SparseArrayProxy;
class SparseArray;

class SparseMatrix
{
    static double               eps;
    RBTree<MatrixIndex, double> tree;
    size_t                      height; 
    size_t                      width;
    size_t                      row;

    public:
    friend SparseArray;
    friend SparseArrayProxy;

    static void set_precision(double);
    static double get_precision();

    SparseMatrix(size_t, size_t);
    SparseMatrix(const SparseMatrix&);
    // ~SparseMatrix(); no need
    SparseMatrix& operator=(const SparseMatrix&);

    double get(size_t, size_t) const;
    void set(size_t, size_t, double);

    size_t num_rows() const;
    size_t num_columns() const;

    SparseArray operator[](size_t);
    SparseArray operator[](size_t) const;
    SparseArray operator*();
    SparseArray operator*() const;
    SparseMatrix& operator+(size_t);
    const SparseMatrix& operator+(size_t) const;
    SparseMatrix& operator-(size_t);
    const SparseMatrix& operator-(size_t) const;

    bool operator==(const SparseMatrix&) const; // after indexation
    bool operator!=(const SparseMatrix&) const; // after indexation

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

    SparseMatrix& operator*=(double val);
    SparseMatrix& operator/=(double val);
    SparseMatrix operator*(double val) const;
    SparseMatrix operator/(double val) const;
};

class SparseArray
{
    bool            modifyable;
    size_t          row;
    size_t          column;
    SparseMatrix    *m;

    Node<MatrixIndex,double> *node;
    
    public:
    friend SparseArrayProxy;
    SparseArray(SparseMatrix *, size_t, bool = true);

    SparseArray& operator+(size_t);
    SparseArray& operator-(size_t);
    SparseArrayProxy operator[](size_t);
    SparseArrayProxy operator*();
};

class SparseArrayProxy
{
    SparseArray   *real;

    public:
    SparseArrayProxy(SparseArray *);
    SparseArrayProxy& operator=(double);
    SparseArrayProxy& operator=(const SparseArrayProxy&);
    
    SparseArrayProxy& perform_operation(void (*)(double&,double), double);
    SparseArrayProxy& operator+=(double);
    SparseArrayProxy& operator-=(double);
    SparseArrayProxy& operator*=(double);
    SparseArrayProxy& operator/=(double);

    bool operator==(double);
    bool operator==(SparseArrayProxy&);
    bool operator!=(double);
    bool operator!=(SparseArrayProxy&);

    operator double() const;
};

#endif