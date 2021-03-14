#ifndef SPARSE_MATRIX_HPP
#define SPARSE_MATRIX_HPP

#include "MatrixIndex.hpp"
#include "RBTree.hpp"
#include <stdexcept>
#include <new>

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

    void check_idx(size_t row) const throw(std::out_of_range);

    friend Cell;
    public:
    static void set_precision(double) throw(std::logic_error);
    static double get_precision() throw();

    SparseMatrix(size_t, size_t) throw(std::logic_error);
    SparseMatrix(const SparseMatrix&) throw(std::bad_alloc);
    ~SparseMatrix() throw();
    SparseMatrix& operator=(const SparseMatrix&) throw(std::bad_alloc);

    /*
    ** Установка значений в ячейках матрицы
    */
    double get(size_t, size_t) const throw(std::out_of_range);
    void set(size_t, size_t, double) throw(std::out_of_range, std::bad_alloc);

    /*
    ** Получение информационных полей
    */
    size_t num_rows() const throw();
    size_t num_columns() const throw();

    /*
    ** Операции разыменования
    */
    ColumnSelector operator[](size_t) throw(std::out_of_range);
    ColumnSelector operator[](size_t) const throw(std::out_of_range);
    ColumnSelector operator*() throw();
    ColumnSelector operator*() const throw();

    /*
    ** Арифметические операции вида:
    ** SparseMatrix.op(SparseMatrix)
    */
    SparseMatrix& perform_operation(void (*)(double&,double), const SparseMatrix&) throw(std::invalid_argument, std::bad_alloc);
    SparseMatrix& operator+=(const SparseMatrix&) throw(std::invalid_argument, std::bad_alloc);
    SparseMatrix& operator-=(const SparseMatrix&) throw(std::invalid_argument, std::bad_alloc);
    SparseMatrix& operator*=(const SparseMatrix&) throw(std::invalid_argument);
    SparseMatrix& operator/=(const SparseMatrix&) throw(std::invalid_argument);
    SparseMatrix dot(const SparseMatrix&) const throw(std::invalid_argument, std::bad_alloc);

    /*
    ** Арифметические операции вида:
    ** SparseMatrix.op(double)
    */
    SparseMatrix& operator*=(double val) throw();
    SparseMatrix& operator/=(double val) throw(std::domain_error);
};

/*
** Логические операции
*/
bool operator==(const SparseMatrix&, const SparseMatrix&) throw();
bool operator!=(const SparseMatrix&, const SparseMatrix&) throw();

/*
** Арифметические операции вида:
** SparseMatrix (*) SparseMatrix
*/
SparseMatrix operator+(const SparseMatrix&, const SparseMatrix&) throw(std::invalid_argument, std::bad_alloc);
SparseMatrix operator-(const SparseMatrix&, const SparseMatrix&) throw(std::invalid_argument, std::bad_alloc);
SparseMatrix operator*(const SparseMatrix&, const SparseMatrix&) throw(std::invalid_argument, std::bad_alloc);
SparseMatrix operator/(const SparseMatrix&, const SparseMatrix&) throw(std::invalid_argument, std::bad_alloc);
SparseMatrix dot(const SparseMatrix&, const SparseMatrix&) throw(std::invalid_argument, std::bad_alloc);

/*
** Арифметичесик операции вида:
** SparseMatrix (*) double или  double (*) SparseMatrix
*/
SparseMatrix operator*(const SparseMatrix&, double) throw(std::bad_alloc);
SparseMatrix operator*(double, const SparseMatrix&) throw(std::bad_alloc);
SparseMatrix operator/(const SparseMatrix&, double) throw(std::domain_error, std::bad_alloc);

/*
** Адресная арифметика
*/
RowSelector operator+(SparseMatrix&, size_t) throw();
RowSelector operator+(const SparseMatrix&, size_t) throw();
RowSelector operator+(size_t, SparseMatrix&) throw();
RowSelector operator+(size_t, const SparseMatrix&) throw();
RowSelector operator-(SparseMatrix&, size_t) throw();
RowSelector operator-(const SparseMatrix&, size_t) throw();

template <class T_C>
class Selector;

template <class T_C>
const Selector<T_C>& operator+(const Selector<T_C>&, size_t) throw();
template <class T_C>
const Selector<T_C>& operator+(size_t, const Selector<T_C>&) throw();
template <class T_C>
const Selector<T_C>& operator-(const Selector<T_C>&, size_t) throw();

template <class T_C>
class Selector
{
    virtual inline void check_idx() const = 0;
protected:
    const SparseMatrix  *matrix;
    const bool          mod;
    size_t              idx;

    explicit Selector(const SparseMatrix *, bool = true, size_t = 0);

    friend const Selector& operator+ <T_C> (const Selector&, size_t);
    friend const Selector& operator+ <T_C> (size_t, const Selector&);
    friend const Selector& operator- <T_C> (const Selector&, size_t);
public:
    virtual T_C     operator[](size_t) const = 0;
    virtual T_C     operator*() const = 0;
};

class RowSelector : private Selector<ColumnSelector>
{
    virtual inline void check_idx() const throw(std::out_of_range);
public:
    explicit RowSelector(const SparseMatrix *, bool = true, size_t = 0) throw();
    virtual ColumnSelector  operator[](size_t) const throw(std::out_of_range);
    virtual ColumnSelector  operator*() const throw(std::out_of_range);
};

class ColumnSelector : public Selector<Cell>
{
    const size_t    row;
    virtual inline void check_idx() const throw(std::out_of_range);
public:
    explicit ColumnSelector(const SparseMatrix *, bool = true, size_t = 0) throw();
    virtual Cell    operator[](size_t) const throw(std::out_of_range);
    virtual Cell    operator*() const throw(std::out_of_range);
};

class Cell
{
    RBTree<MatrixIndex,double>&         tree;
    const MatrixIndex                   position;
    const bool                          mod;
    Node<MatrixIndex,double> *          node;

public:
    Cell(const SparseMatrix *, size_t, size_t, bool) throw();

    Cell& operator=(double) throw(std::logic_error);
    Cell& operator=(const Cell&) throw(std::logic_error);
    
    Cell& perform_operation(void (*)(double&,double), double) throw(std::logic_error, std::bad_alloc);
    Cell& operator+=(double) throw(std::logic_error, std::bad_alloc);
    Cell& operator-=(double) throw(std::logic_error, std::bad_alloc);
    Cell& operator*=(double) throw(std::logic_error);
    Cell& operator/=(double) throw(std::logic_error, std::domain_error);

    bool operator==(double) const throw();
    bool operator==(Cell&) const throw();
    bool operator!=(double) const throw();
    bool operator!=(Cell&) const throw();

    operator double() const;
};

#endif