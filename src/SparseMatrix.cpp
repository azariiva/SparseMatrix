#include "SparseMatrix.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

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
    return ColumnSelector(RowSelector(this, true, row));
}

ColumnSelector SparseMatrix::operator[](size_t row) const {
    return ColumnSelector(RowSelector(const_cast<SparseMatrix *>(this), false, row), false);
}

ColumnSelector SparseMatrix::operator*() {
   return ColumnSelector(RowSelector(this, true));
}

ColumnSelector SparseMatrix::operator*() const {
    return ColumnSelector(RowSelector(const_cast<SparseMatrix *>(this), false), false);
}

RowSelector SparseMatrix::operator+(size_t row) {
    return RowSelector(this, true, row);
}

RowSelector SparseMatrix::operator+(size_t row) const {
    return RowSelector(const_cast<SparseMatrix *>(this), false, row);
}

RowSelector SparseMatrix::operator-(size_t row) {
    return RowSelector(this, true, -row);
}

RowSelector SparseMatrix::operator-(size_t row) const {
    return RowSelector(const_cast<SparseMatrix *>(this), -row);
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

SparseMatrix SparseMatrix::operator+(const SparseMatrix& rv) const {
    SparseMatrix m(*this);
    m += rv;
    return m;
}

SparseMatrix SparseMatrix::operator-(const SparseMatrix& rv) const {
    SparseMatrix m(*this);
    m -= rv;
    return m;
}

SparseMatrix SparseMatrix::operator*(const SparseMatrix& rv) const {
    SparseMatrix m(*this);
    m *= rv;
    return m;
}

SparseMatrix SparseMatrix::operator/(const SparseMatrix& rv) const {
    SparseMatrix m(*this);
    m /= rv;
    return m;
}

SparseMatrix SparseMatrix::dot(const SparseMatrix& rv) const {
    // TODO: rewrite after transposed matrix is implemented
    SparseMatrix m(height, rv.width);

    if (width != rv.height) {
        throw std::invalid_argument("SparseMatrix::dot : operands could not be broadcast together");
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

SparseMatrix SparseMatrix::operator*(double val) const {
    SparseMatrix m(*this);
    m *= val;
    return m;
}

SparseMatrix SparseMatrix::operator/(double val) const {
    return operator*(1.0 / val);
}

bool SparseMatrix::operator==(const SparseMatrix& rv) const {
    if (height != rv.height || width != rv.width) {
        return false;
    }
    for (size_t row_ = 0; row_ < height; row_++) {
        for (size_t column = 0; column < width; column++) {
            if (operator[](row_)[column] != rv[row_][column]) {
                return false;
            }
        }
    }
    return true;
}

bool SparseMatrix::operator!=(const SparseMatrix& rv) const {
    if (height != rv.height || width != rv.width) {
        return true;
    }
    for (size_t row_ = 0; row_ < height; row_++) {
        for (size_t column = 0; column < width; column++) {
            if (operator[](row_)[column] == rv[row_][column]) {
                return false;
            }
        }
    }
    return true;
}

RowSelector operator+(size_t lv, SparseMatrix& rv) {
    return rv.operator+(lv);
}

RowSelector operator+(size_t lv, const SparseMatrix& rv) {
    return rv.operator+(lv);
}

SparseMatrix operator*(double lv, const SparseMatrix& rv) {
    return rv.operator*(lv);
}

double SparseMatrix::eps  = 1e-8;
size_t SparseMatrix::instance_quantity = 0;

template <class T_P, class T_C>
Selector<T_P,T_C>::Selector(T_P p_, bool mod_, size_t idx_) : p(p_) {
    mod = mod_;
    idx = idx_;
}

template <class T_P, class T_C>
const Selector<T_P,T_C>& Selector<T_P,T_C>::operator+(size_t i) const {
    const_cast<Selector<T_P,T_C>*>(this)->idx += i;
    return *this;
}

template <class T_P, class T_C>
const Selector<T_P,T_C>& Selector<T_P,T_C>::operator-(size_t i) const {
    const_cast<Selector<T_P,T_C>*>(this)->idx -= i;
    return *this;
}

template <class T_P, class T_C>
const Selector<T_P,T_C>& operator+(size_t i, const Selector<T_P,T_C>& s) {
    return s.operator+(i);
}

RowSelector::RowSelector(SparseMatrix* p_, bool mod_, size_t idx_) : Selector<SparseMatrix *,ColumnSelector>(p_, mod_, idx_) {}

void RowSelector::check_idx() const {
    if (idx > p->num_columns()) {
        //TODO: rewrite error message
        throw std::out_of_range("out of range");
    }
}

ColumnSelector RowSelector::operator[](size_t i) const {
    const_cast<RowSelector *>(this)->idx += i;
    check_idx();
    return ColumnSelector(*this, mod);
}

ColumnSelector RowSelector::operator*() const {
    check_idx();
    return ColumnSelector(*this, mod);
}

ColumnSelector::ColumnSelector(RowSelector p_, bool mod_, size_t idx_) : Selector<RowSelector,Cell>(p_, mod_, idx_) {}

void ColumnSelector::check_idx() const {
    if (idx > p.p->width) {
        //TODO: rewrite error message
        throw std::out_of_range("out of range");
    }
}

Cell ColumnSelector::operator[](size_t i) const {
    const_cast<ColumnSelector *>(this)->idx += i;
    check_idx();
    return Cell(p.p->tree, p.idx, idx, mod);
}

Cell ColumnSelector::operator*() const {
    check_idx();
    return Cell(p.p->tree, p.idx, idx, mod);
}

Cell::Cell(RBTree<MatrixIndex,double>& tree_, size_t row, size_t column, bool mod_) : position(row, column), mod(mod_), tree(tree_) {
    node = tree.get_node(position);
}

Cell& Cell::operator=(double rv) {
    if (mod == false) {
        // TODO: rewrite error message
        throw std::logic_error("could not modify constant value");
    }
    if (node == Node<MatrixIndex,double>::nil_node) {
        if (fabs(rv - 0.0) > SparseMatrix::get_precision()) {
            node = tree.insert(position, rv);
        }
    } else {
        if (fabs(rv - 0.0) > SparseMatrix::get_precision()) {
            node->item = rv;
        } else {
            tree.remove(node);
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
            node = tree.insert(position, tmp);
        }
    } else if (node_exist) {
        tree.remove(node);
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

const RowSelector& operator+(size_t lv, RowSelector& rv) {
    const_cast<RowSelector*>(&rv)->idx += lv;
    return rv;
}

const ColumnSelector& operator+(size_t lv, const ColumnSelector& rv) {
    const_cast<ColumnSelector*>(&rv)->idx += lv;
    return rv;
}

template class Selector<RowSelector, Cell>;
template class Selector<SparseMatrix, ColumnSelector>;