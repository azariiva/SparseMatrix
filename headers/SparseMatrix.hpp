#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "MatrixIndex.hpp"
#include "RBTree.hpp"

class SparseMatrixRowProxy;
class SparseMatrixRow;

class RowSelector;
class ColumnSelector;

class SparseMatrix
{
    static double               eps;
    static size_t               instance_quantity;

    RBTree<MatrixIndex, double> tree;
    size_t                      height; 
    size_t                      width;

    friend RowSelector operator+(size_t, SparseMatrix&);
    friend RowSelector operator+(size_t, const SparseMatrix&);
    friend SparseMatrix operator*(double, const SparseMatrix&);
    
    friend ColumnSelector;
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

    ColumnSelector operator[](size_t);
    ColumnSelector operator[](size_t) const;
    ColumnSelector operator*();
    ColumnSelector operator*() const;
    RowSelector operator+(size_t);
    RowSelector operator+(size_t) const;
    RowSelector operator-(size_t);
    RowSelector operator-(size_t) const;

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

template <class T_P, class T_C>
class Selector
{
    virtual inline void check_idx() const = 0;
protected:
    T_P      p;
    bool    mod;
    size_t  idx;

    Selector(T_P, bool = true, size_t = 0);

    friend const RowSelector& operator+(size_t lv, RowSelector& rv);
    friend const ColumnSelector& operator+(size_t lv, const ColumnSelector& rv);
public:
    const Selector& operator+(size_t) const;
    const Selector& operator-(size_t) const;
    virtual T_C     operator[](size_t) const = 0;
    virtual T_C     operator*() const = 0;
};

class ColumnSelector;

class RowSelector : public Selector<SparseMatrix *, ColumnSelector>
{
    virtual inline void check_idx() const;

    friend ColumnSelector;
public:

    RowSelector(SparseMatrix *, bool = true, size_t = 0);
    virtual ColumnSelector  operator[](size_t) const;
    virtual ColumnSelector  operator*() const;
};

class Cell;

class ColumnSelector : public Selector<RowSelector, Cell>
{
    virtual inline void check_idx() const;
public:
    ColumnSelector(RowSelector, bool = true, size_t = 0);
    virtual Cell    operator[](size_t) const;
    virtual Cell    operator*() const;
};

class Cell
{
    Node<MatrixIndex,double>    *node;
    MatrixIndex                 position;
    bool                        mod;
    RBTree<MatrixIndex,double>  &tree;

public:
    Cell(RBTree<MatrixIndex,double>&, size_t, size_t, bool);

    Cell& operator=(double);
    Cell& operator=(const Cell&);
    
    Cell& perform_operation(void (*)(double&,double), double);
    Cell& operator+=(double);
    Cell& operator-=(double);
    Cell& operator*=(double);
    Cell& operator/=(double);

    bool operator==(double) const;
    bool operator==(Cell&) const;
    bool operator!=(double) const;
    bool operator!=(Cell&) const;

    operator double() const;
};

#endif