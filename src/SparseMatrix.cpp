#include "SparseMatrix.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

double SparseMatrix::eps  = 1e-8;
size_t SparseMatrix::instance_quantity = 0;

void SparseMatrix::set_precision(double eps_) {
    if (SparseMatrix::instance_quantity != 0) {
        throw std::runtime_error("SparseMatrix::set_precision : precision could not be modified when instances of SparseMatrix exist");
    } else {
        SparseMatrix::eps = eps_;
    }
}

double SparseMatrix::get_precision() {
    return SparseMatrix::eps;
}

SparseMatrix::SparseMatrix(size_t height_, size_t width_) {
    height = height_;
    width = width_;
    SparseMatrix::instance_quantity++;
}

SparseMatrix::SparseMatrix(const SparseMatrix& src) : tree(src.tree) {
    height = src.height;
    width = src.width;
    SparseMatrix::instance_quantity++;
}

SparseMatrix::~SparseMatrix() {
    SparseMatrix::instance_quantity--;
}

SparseMatrix& SparseMatrix::operator=(const SparseMatrix& src) {
    if (this == &src) {
        return *this;
    }
    height = src.height;
    width = src.width;
    tree = src.tree;
    return *this;
}

double SparseMatrix::get(size_t row, size_t column) const {
    Node<MatrixIndex, double> *node;

    if (row >= height || column >= width) {
        throw std::out_of_range("SparseMatrix::get : out of bounds exception(height or width exceeded)");
    }
    node = tree.get_node(MatrixIndex(row, column));
    if (node == RBTree<MatrixIndex, double>::nil_node) {
        return 0.0;
    }
    return node->item;
}

void SparseMatrix::set(size_t row, size_t column, double val) {
    if (row >= height || column >= width) {
        throw std::out_of_range("SparseMatrix::set : out of bounds exception(height or width exceeded)");
    }
    if (fabs(val - 0.0) <= SparseMatrix::eps) {
        tree.remove(MatrixIndex(row, column));
    } else {
        tree.insert(MatrixIndex(row, column), val);
    }
}

size_t SparseMatrix::num_rows() const {
    return height;
}

size_t SparseMatrix::num_columns() const {
    return width;
}

ColumnSelector SparseMatrix::operator[](size_t row) {
    return ColumnSelector(this, row);
}

ColumnSelector SparseMatrix::operator[](size_t row) const {
    return ColumnSelector(this, row, false);
}

ColumnSelector SparseMatrix::operator*() {
   return ColumnSelector(this);
}

ColumnSelector SparseMatrix::operator*() const {
    return ColumnSelector(this, 0, false);
}

SparseMatrix& SparseMatrix::perform_operation(void (*op)(double&,double), const SparseMatrix&rv) {
    // TODO: optimize a little
    if (height != rv.height || width != rv.width) {
        throw std::invalid_argument("SparseMatrix::perform_operation : operands could not be broadcast together");
    }
    for (size_t row = 0; row < height; row++) {
        for (size_t column = 0; column < width; column++) {
            (operator[](row)[column]).perform_operation(op, rv[row][column]);
        }
    }
    return *this;
}

/*
    Все 4 нижестоящие функции - обёртки вокруг SparseMatrixRowProxy::perform_operation()
*/
SparseMatrix& SparseMatrix::operator+=(const SparseMatrix& rv) {
    return perform_operation([](double &lv, double rv_) {lv += rv_;}, rv);
}

SparseMatrix& SparseMatrix::operator-=(const SparseMatrix& rv) {
    return perform_operation([](double &lv, double rv_) {lv -= rv_;}, rv);
}

SparseMatrix& SparseMatrix::operator*=(const SparseMatrix& rv) {
    return perform_operation([](double &lv, double rv_) {lv *= rv_;}, rv);
}

SparseMatrix& SparseMatrix::operator/=(const SparseMatrix& rv) {
    return perform_operation([](double &lv, double rv_) {lv /= rv_;}, rv);
}

SparseMatrix& SparseMatrix::operator*=(double val) {
    for (size_t row_ = 0; row_ < height; row_++) {
        for (size_t column = 0; column < height; column++) {
            operator[](row_)[column] *= val;
        }
    }
    return *this;
}

SparseMatrix& SparseMatrix::operator/=(double val) {
    return operator*=(1.0 / val);
}

SparseMatrix SparseMatrix::dot(const SparseMatrix& rv) const {
    // TODO: rewrite after transposed matrix is implemented
    SparseMatrix m(height, rv.width);

    if (width != rv.height) {
        throw std::invalid_argument("dot : operands could not be broadcast together");
    }
    for (size_t row_ = 0; row_ < height; row_++) {
        for (size_t column = 0; column < rv.width; column++) {
            for (size_t k = 0; k < width; k++) {
                m[row_][column] += operator[](row_)[k] * rv[k][column];
            }
        }
    }
    return m;
}

bool operator==(const SparseMatrix& lv, const SparseMatrix& rv) {
    if (lv.num_rows() != rv.num_rows() || lv.num_columns() != rv.num_columns()) {
        return false;
    }
    for (size_t row = 0; row < lv.num_rows(); row++) {
        for (size_t column = 0; column < lv.num_columns(); column++) {
            if (lv.get(row, column) != rv.get(row, column)) {
                return false;
            }
        }
    }
    return true;
}

bool operator!=(const SparseMatrix& lv, const SparseMatrix& rv) {
    return !(lv == rv);
}

SparseMatrix operator+(const SparseMatrix& lv, const SparseMatrix& rv) {
    SparseMatrix res(lv);
    res += rv;
    return res;
}

SparseMatrix operator-(const SparseMatrix& lv, const SparseMatrix& rv) {
    SparseMatrix res(lv);
    res -= rv;
    return res;
}

SparseMatrix operator*(const SparseMatrix& lv, const SparseMatrix& rv) {
    SparseMatrix res(lv);
    res *= rv;
    return res;
}

SparseMatrix operator/(const SparseMatrix& lv, const SparseMatrix& rv) {
    SparseMatrix res(lv);
    res /= rv;
    return res;
}

SparseMatrix dot(const SparseMatrix& lv, const SparseMatrix& rv) {
    return lv.dot(rv);
}

SparseMatrix operator*(const SparseMatrix& lv, double rv) {
    if (fabs(rv - 0.0) < SparseMatrix::get_precision()) {
        return SparseMatrix(lv.num_rows(), lv.num_columns());
    }
    SparseMatrix result(lv);
    result *= rv;
    return result;
}

SparseMatrix operator*(double lv, const SparseMatrix& rv) {
    return operator*(rv, lv);
}

SparseMatrix operator/(const SparseMatrix& lv, double rv) {
    return operator*(lv, 1.0 / rv);
}

RowSelector operator+(SparseMatrix& matrix, size_t row) {
    return RowSelector(&matrix, true, row);
}

RowSelector operator+(const SparseMatrix& matrix, size_t row) {
    return RowSelector(&matrix, false, row);
}

RowSelector operator+(size_t row, SparseMatrix& matrix) {
    return RowSelector(&matrix, true, row);
}

RowSelector operator+(size_t row, const SparseMatrix& matrix) {
    return RowSelector(&matrix, false, row);
}

RowSelector operator-(SparseMatrix& matrix, size_t row) {
    return RowSelector(&matrix, true, -row);
}

RowSelector operator-(const SparseMatrix& matrix, size_t row) {
    return RowSelector(&matrix, false, -row);
}

template <class T_C>
Selector<T_C>::Selector(const SparseMatrix *matrix_, bool mod_, size_t idx_) : matrix(matrix_), mod(mod_), idx(idx_) {}

template <class T_C>
const Selector<T_C>& operator+(const Selector<T_C>& lv, size_t rv) {
    const_cast<Selector<T_C> *>(&lv)->idx += rv;
    return lv;
}

template <class T_C>
const Selector<T_C>& operator+(size_t lv, const Selector<T_C>& rv) {
    const_cast<Selector<T_C> *>(&rv)->idx += lv;
    return rv;
}

template <class T_C>
const Selector<T_C>& operator-(const Selector<T_C>& lv, size_t rv) {
    const_cast<Selector<T_C> *>(&lv)->idx -= rv;
    return lv;
}

RowSelector::RowSelector(const SparseMatrix* matrix_, bool mod_, size_t idx_) : Selector<ColumnSelector>(matrix_, mod_, idx_) {}

void RowSelector::check_idx() const {
    if (idx > matrix->num_columns()) {
        //TODO: rewrite error message
        throw std::out_of_range("out of range");
    }
}

ColumnSelector RowSelector::operator[](size_t i) const {
    const_cast<RowSelector *>(this)->idx += i;
    check_idx();
    return ColumnSelector(matrix, idx, mod);
}

ColumnSelector RowSelector::operator*() const {
    check_idx();
    return ColumnSelector(matrix, idx, mod);
}

ColumnSelector::ColumnSelector(const SparseMatrix *matrix_, size_t row_, bool mod_, size_t idx_) : Selector<Cell>(matrix_, mod_, idx_), row(row_) {}

void ColumnSelector::check_idx() const {
    if (idx > matrix->num_columns()) {
        //TODO: rewrite error message
        throw std::out_of_range("out of range");
    }
}

Cell ColumnSelector::operator[](size_t i) const {
    const_cast<ColumnSelector *>(this)->idx += i;
    check_idx();
    return Cell(matrix, row, idx, mod);
}

Cell ColumnSelector::operator*() const {
    check_idx();
    return Cell(matrix, row, idx, mod);
}

Cell::Cell(const SparseMatrix *matrix_, size_t row, size_t column, bool mod_) : matrix(const_cast<SparseMatrix *>(matrix_)), position(row, column), mod(mod_) {
    node = matrix->tree.get_node(position);
}

Cell& Cell::operator=(double rv) {
    if (mod == false) {
        // TODO: rewrite error message
        throw std::logic_error("could not modify constant value");
    }
    if (node == Node<MatrixIndex,double>::nil_node) {
        if (fabs(rv - 0.0) > SparseMatrix::get_precision()) {
            node = matrix->tree.insert(position, rv);
        }
    } else {
        if (fabs(rv - 0.0) > SparseMatrix::get_precision()) {
            node->item = rv;
        } else {
            matrix->tree.remove(node);
            node = Node<MatrixIndex,double>::nil_node;
        }
    }
    return *this;
}

Cell::operator double() const {
    if (node == Node<MatrixIndex,double>::nil_node) {
        return 0.0;
    }
    return node->item;
}

Cell& Cell::operator=(const Cell& rv) {
    if (this == &rv) {
        return *this;
    }
    return operator=(double(rv));
}

Cell& Cell::perform_operation(void (*op)(double&,double), double rv) {
    bool node_exist = (node != Node<MatrixIndex,double>::nil_node);
    double tmp = (node_exist ? node->item : 0.0);

    if (mod == false) {
        // TODO: change error message
        throw std::runtime_error("SparseMatrixRowProxy::perform_operation : constant instance could not be modified");
    }
    op(tmp, rv);
    if (fabs(tmp - 0.0) > SparseMatrix::get_precision()) {
        if (node_exist) {
            node->item = tmp;
        } else {
            node = matrix->tree.insert(position, tmp);
        }
    } else if (node_exist) {
        matrix->tree.remove(node);
        node = Node<MatrixIndex,double>::nil_node;
    }
    return *this;
}

Cell& Cell::operator+=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv += rv_;}, rv);
}

Cell& Cell::operator-=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv -= rv_;}, rv);
}

Cell& Cell::operator*=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv *= rv_;}, rv);
}

Cell& Cell::operator/=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv /= rv_;}, rv);
}

bool Cell::operator==(double rv) const {
    return (fabs(operator double() - rv) < SparseMatrix::get_precision());
}

bool Cell::operator==(Cell& rv) const {
    return (fabs(operator double() - double(rv)) < SparseMatrix::get_precision());
}

bool Cell::operator!=(double rv) const {
    return (fabs(operator double() - rv) >= SparseMatrix::get_precision());
}

bool Cell::operator!=(Cell& rv) const {
    return (fabs(operator double() - double(rv)) >= SparseMatrix::get_precision());
}

template class Selector<Cell>;
template class Selector<ColumnSelector>;
template const Selector<ColumnSelector>& operator+(const Selector<ColumnSelector>&, size_t);
template const Selector<ColumnSelector>& operator+(size_t, const Selector<ColumnSelector>&);
template const Selector<ColumnSelector>& operator-(const Selector<ColumnSelector>&, size_t);

template const Selector<Cell>& operator+(const Selector<Cell>&, size_t);
template const Selector<Cell>& operator+(size_t, const Selector<Cell>&);
template const Selector<Cell>& operator-(const Selector<Cell>&, size_t);

