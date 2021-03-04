#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "MatrixIndex.hpp"
#include "RBTree.hpp"

class Cell;
class ColumnSelector;
class RowSelector;

class SparseMatrix
{
    static double               eps;
    static size_t               instance_quantity;

    RBTree<MatrixIndex, double> tree;
    size_t                      height; 
    size_t                      width;

    void check_idx(size_t row) const;
    
    friend Cell;
    public:
    static void set_precision(double);
    static double get_precision();

    SparseMatrix(size_t, size_t);
    SparseMatrix(const SparseMatrix&);
    ~SparseMatrix();
    SparseMatrix& operator=(const SparseMatrix&);

    /*
    ** Установка значений в ячейках матрицы
    */
    double get(size_t, size_t) const;
    void set(size_t, size_t, double);

    /*
    ** Получение информационных полей
    */
    size_t num_rows() const;
    size_t num_columns() const;

    /*
    ** Операции разыменования
    */
    ColumnSelector operator[](size_t);
    ColumnSelector operator[](size_t) const;
    ColumnSelector operator*();
    ColumnSelector operator*() const;

    /*
    ** Арифметические операции вида:
    ** SparseMatrix.op(SparseMatrix)
    */
    SparseMatrix& perform_operation(void (*)(double&,double), const SparseMatrix&);
    SparseMatrix& operator+=(const SparseMatrix&);
    SparseMatrix& operator-=(const SparseMatrix&);
    SparseMatrix& operator*=(const SparseMatrix&);
    SparseMatrix& operator/=(const SparseMatrix&);
    SparseMatrix dot(const SparseMatrix&) const;

    /*
    ** Арифметические операции вида:
    ** SparseMatrix.op(double)
    */
    SparseMatrix& operator*=(double val);
    SparseMatrix& operator/=(double val);
};

/*
** Логические операции
*/
bool operator==(const SparseMatrix&, const SparseMatrix&);
bool operator!=(const SparseMatrix&, const SparseMatrix&);

/*
** Арифметические операции вида:
** SparseMatrix (*) SparseMatrix
*/
SparseMatrix operator+(const SparseMatrix&, const SparseMatrix&);
SparseMatrix operator-(const SparseMatrix&, const SparseMatrix&);
SparseMatrix operator*(const SparseMatrix&, const SparseMatrix&);
SparseMatrix operator/(const SparseMatrix&, const SparseMatrix&);
SparseMatrix dot(const SparseMatrix&, const SparseMatrix&);

/*
** Арифметичесик операции вида:
** SparseMatrix (*) double или  double (*) SparseMatrix
*/
SparseMatrix operator*(const SparseMatrix&, double);
SparseMatrix operator*(double, const SparseMatrix&);
SparseMatrix operator/(const SparseMatrix&, double);

/*
** Адресная арифметика
*/
RowSelector operator+(SparseMatrix&, size_t);
RowSelector operator+(const SparseMatrix&, size_t);
RowSelector operator+(size_t, SparseMatrix&);
RowSelector operator+(size_t, const SparseMatrix&);
RowSelector operator-(SparseMatrix&, size_t);
RowSelector operator-(const SparseMatrix&, size_t);

template <class T_C>
class Selector;

template <class T_C>
const Selector<T_C>& operator+(const Selector<T_C>&, size_t);
template <class T_C>
const Selector<T_C>& operator+(size_t, const Selector<T_C>&);
template <class T_C>
const Selector<T_C>& operator-(const Selector<T_C>&, size_t);

template <class T_C>
class Selector
{
    virtual inline void check_idx() const = 0;
protected:
    const SparseMatrix  *matrix;
    const bool          mod;
    size_t              idx;

    Selector(const SparseMatrix *, bool = true, size_t = 0);

    friend const Selector& operator+ <T_C> (const Selector&, size_t);
    friend const Selector& operator+ <T_C> (size_t, const Selector&);
    friend const Selector& operator- <T_C> (const Selector&, size_t);
public:
    virtual T_C     operator[](size_t) const = 0;
    virtual T_C     operator*() const = 0;
};

class RowSelector : public Selector<ColumnSelector>
{
    virtual inline void check_idx() const;

    friend ColumnSelector;
public:
    explicit RowSelector(const SparseMatrix *, bool = true, size_t = 0);
    virtual ColumnSelector  operator[](size_t) const;
    virtual ColumnSelector  operator*() const;
};

// const RowSelector& operator+(const RowSelector&, size_t);
// const RowSelector& operator+(size_t, const RowSelector&);
// const RowSelector& operator-(const RowSelector&, size_t);

class ColumnSelector : public Selector<Cell>
{
    const size_t    row;
    virtual inline void check_idx() const;
public:
    explicit ColumnSelector(const SparseMatrix *, size_t = 0, bool = true, size_t = 0);
    virtual Cell    operator[](size_t) const;
    virtual Cell    operator*() const;
};

// const ColumnSelector& operator+(const ColumnSelector&, size_t);
// const ColumnSelector& operator+(size_t, const ColumnSelector&);
// const ColumnSelector& operator-(const ColumnSelector&, size_t);

class Cell
{
    SparseMatrix * const                matrix;
    const MatrixIndex                   position;
    const bool                          mod;
    Node<MatrixIndex,double> *          node;

public:
    Cell(const SparseMatrix *, size_t, size_t, bool);

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