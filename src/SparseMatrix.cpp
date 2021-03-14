#include "SparseMatrix.hpp"
#include <cmath>
#include <iostream>

double SparseMatrix::eps  = 1e-8;
size_t SparseMatrix::instance_quantity = 0;

void SparseMatrix::check_idx(size_t row) const throw(std::out_of_range) {
    if (row >= height) {
        throw std::out_of_range("out of bounds exception (height exceeded)");
    }
}

void SparseMatrix::set_precision(double eps_) throw(std::logic_error) {
    if (SparseMatrix::instance_quantity != 0) {
        throw std::logic_error("precision could not be modified when instances exist");
    } else {
        eps = eps_;
    }
}

double SparseMatrix::get_precision() throw() {
    return eps;
}

SparseMatrix::SparseMatrix(size_t height_, size_t width_) throw(std::logic_error) : height(height_), width(width_) {
    if (height == 0 || width == 0) {
        throw std::logic_error("height and width of matrix could not be 0");
    }
    SparseMatrix::instance_quantity++;
}

SparseMatrix::SparseMatrix(const SparseMatrix& src) throw(std::bad_alloc) : tree(src.tree), height(src.height), width(src.width) {
    SparseMatrix::instance_quantity++;
}

SparseMatrix::~SparseMatrix() throw() {
    SparseMatrix::instance_quantity--;
}

SparseMatrix& SparseMatrix::operator=(const SparseMatrix& src) throw(std::bad_alloc) {
    if (this == &src) {
        return *this;
    }
    height = src.height;
    width = src.width;
    tree = src.tree;
    return *this;
}

double SparseMatrix::get(size_t row, size_t column) const throw(std::out_of_range) {
    Node<MatrixIndex, double> *node;

    if (row >= height || column >= width) {
        throw std::out_of_range("out of bounds exception (height or width exceeded)");
    }
    node = tree.get_node(MatrixIndex(row, column));
    if (node == RBTree<MatrixIndex, double>::nil_node) {
        return 0.0;
    }
    return node->item;
}

void SparseMatrix::set(size_t row, size_t column, double val) throw(std::out_of_range, std::bad_alloc) {
    if (row >= height || column >= width) {
        throw std::out_of_range("out of bounds exception (height or width exceeded)");
    }
    if (fabs(val - 0.0) <= eps) {
        tree.remove(MatrixIndex(row, column));
    } else {
        tree.insert(MatrixIndex(row, column), val);
    }
}

size_t SparseMatrix::num_rows() const throw() {
    return height;
}

size_t SparseMatrix::num_columns() const throw() {
    return width;
}

ColumnSelector SparseMatrix::operator[](size_t row) throw(std::out_of_range) {
    check_idx(row);
    return ColumnSelector(this, true, row);
}

ColumnSelector SparseMatrix::operator[](size_t row) const throw(std::out_of_range) {
    check_idx(row);
    return ColumnSelector(this, false, row);
}

ColumnSelector SparseMatrix::operator*() throw() {
   return ColumnSelector(this);
}

ColumnSelector SparseMatrix::operator*() const throw() {
    return ColumnSelector(this, false);
}

SparseMatrix& SparseMatrix::perform_operation(void (*op)(double&,double), const SparseMatrix&rv) throw(std::invalid_argument, std::bad_alloc) {
    // TODO: optimize a little
    if (height != rv.height || width != rv.width) {
        throw std::invalid_argument("operands could not be broadcast together");
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
SparseMatrix& SparseMatrix::operator+=(const SparseMatrix& rv) throw(std::invalid_argument, std::bad_alloc) {
    return perform_operation([](double &lv, double rv_) {lv += rv_;}, rv);
}

SparseMatrix& SparseMatrix::operator-=(const SparseMatrix& rv) throw(std::invalid_argument, std::bad_alloc) {
    return perform_operation([](double &lv, double rv_) {lv -= rv_;}, rv);
}

SparseMatrix& SparseMatrix::operator*=(const SparseMatrix& rv) throw(std::invalid_argument) {
    return perform_operation([](double &lv, double rv_) {lv *= rv_;}, rv);
}

SparseMatrix& SparseMatrix::operator/=(const SparseMatrix& rv) throw(std::invalid_argument) {
    return perform_operation([](double &lv, double rv_) {lv /= rv_;}, rv);
}

SparseMatrix& SparseMatrix::operator*=(double val) throw() {
    //TODO: clean all tree if val <= eps
    for (size_t row_ = 0; row_ < height; row_++) {
        for (size_t column = 0; column < height; column++) {
            operator[](row_)[column] *= val;
        }
    }
    return *this;
}

SparseMatrix& SparseMatrix::operator/=(double val) throw(std::domain_error) {
    if (fabs(val - 0.0) <= eps) {
        throw std::domain_error("could not divide by 0");
    }
    return operator*=(1.0 / val);
}

SparseMatrix SparseMatrix::dot(const SparseMatrix& rv) const throw(std::invalid_argument, std::bad_alloc) {
    // TODO: rewrite after transposed matrix is implemented
    SparseMatrix m(height, rv.width);

    if (width != rv.height) {
        throw std::invalid_argument("operands could not be broadcast together");
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

bool operator==(const SparseMatrix& lv, const SparseMatrix& rv) throw() {
    if (lv.num_rows() != rv.num_rows() || lv.num_columns() != rv.num_columns()) {
        return false;
    }
    for (size_t row = 0; row < lv.num_rows(); row++) {
        for (size_t column = 0; column < lv.num_columns(); column++) {
            if (fabs(lv.get(row, column) - rv.get(row, column)) > SparseMatrix::get_precision()) {
                return false;
            }
        }
    }
    return true;
}

bool operator!=(const SparseMatrix& lv, const SparseMatrix& rv) throw() {
    return !(lv == rv);
}

SparseMatrix operator+(const SparseMatrix& lv, const SparseMatrix& rv) throw(std::invalid_argument, std::bad_alloc) {
    SparseMatrix res(lv);
    res += rv;
    return res;
}

SparseMatrix operator-(const SparseMatrix& lv, const SparseMatrix& rv) throw(std::invalid_argument, std::bad_alloc) {
    SparseMatrix res(lv);
    res -= rv;
    return res;
}

SparseMatrix operator*(const SparseMatrix& lv, const SparseMatrix& rv) throw(std::invalid_argument, std::bad_alloc) {
    SparseMatrix res(lv);
    res *= rv;
    return res;
}

SparseMatrix operator/(const SparseMatrix& lv, const SparseMatrix& rv) throw(std::invalid_argument, std::bad_alloc) {
    SparseMatrix res(lv);
    res /= rv;
    return res;
}

SparseMatrix dot(const SparseMatrix& lv, const SparseMatrix& rv) throw(std::invalid_argument, std::bad_alloc) {
    return lv.dot(rv);
}

SparseMatrix operator*(const SparseMatrix& lv, double rv) throw(std::bad_alloc) {
    if (fabs(rv - 0.0) <= SparseMatrix::get_precision()) {
        return SparseMatrix(lv.num_rows(), lv.num_columns());
    }
    SparseMatrix result(lv);
    result *= rv;
    return result;
}

SparseMatrix operator*(double lv, const SparseMatrix& rv) throw(std::bad_alloc) {
    return operator*(rv, lv);
}

SparseMatrix operator/(const SparseMatrix& lv, double rv) throw(std::domain_error, std::bad_alloc) {
    if (fabs(rv - 0.0) <= SparseMatrix::get_precision()) {
        throw std::domain_error("could not divide by 0");
    }
    return operator*(lv, 1.0 / rv);
}

RowSelector operator+(SparseMatrix& matrix, size_t row) throw() {
    return RowSelector(&matrix, true, row);
}

RowSelector operator+(const SparseMatrix& matrix, size_t row) throw() {
    return RowSelector(&matrix, false, row);
}

RowSelector operator+(size_t row, SparseMatrix& matrix) throw() {
    return RowSelector(&matrix, true, row);
}

RowSelector operator+(size_t row, const SparseMatrix& matrix) throw() {
    return RowSelector(&matrix, false, row);
}

RowSelector operator-(SparseMatrix& matrix, size_t row) throw() {
    return RowSelector(&matrix, true, -row);
}

RowSelector operator-(const SparseMatrix& matrix, size_t row) throw() {
    return RowSelector(&matrix, false, -row);
}

template <class T_C>
Selector<T_C>::Selector(const SparseMatrix *matrix_, bool mod_, size_t idx_) : matrix(matrix_), mod(mod_), idx(idx_) {}

template <class T_C>
const Selector<T_C>& operator+(const Selector<T_C>& lv, size_t rv) throw() {
    const_cast<Selector<T_C> *>(&lv)->idx += rv;
    return lv;
}

template <class T_C>
const Selector<T_C>& operator+(size_t lv, const Selector<T_C>& rv) throw() {
    const_cast<Selector<T_C> *>(&rv)->idx += lv;
    return rv;
}

template <class T_C>
const Selector<T_C>& operator-(const Selector<T_C>& lv, size_t rv) throw() {
    const_cast<Selector<T_C> *>(&lv)->idx -= rv;
    return lv;
}

RowSelector::RowSelector(const SparseMatrix* matrix_, bool mod_, size_t idx_) throw() : Selector<ColumnSelector>(matrix_, mod_, idx_) {}

void RowSelector::check_idx() const throw(std::out_of_range) {
    if (idx >= matrix->num_columns()) {
        throw std::out_of_range("out of bounds exception(height exceeded)");
    }
}

ColumnSelector RowSelector::operator[](size_t i) const throw(std::out_of_range) {
    const_cast<RowSelector *>(this)->idx += i;
    check_idx();
    return ColumnSelector(matrix, mod, idx);
}

ColumnSelector RowSelector::operator*() const throw(std::out_of_range) {
    check_idx();
    return ColumnSelector(matrix, mod, idx);
}

ColumnSelector::ColumnSelector(const SparseMatrix *matrix_, bool mod_, size_t row_) throw() : Selector<Cell>(matrix_, mod_), row(row_) {}

void ColumnSelector::check_idx() const throw(std::out_of_range) {
    if (idx >= matrix->num_columns()) {
        throw std::out_of_range("out of bounds exception(width exceeded)");
    }
}

Cell ColumnSelector::operator[](size_t i) const throw(std::out_of_range) {
    const_cast<ColumnSelector *>(this)->idx += i;
    check_idx();
    return Cell(matrix, row, idx, mod);
}

Cell ColumnSelector::operator*() const throw(std::out_of_range) {
    check_idx();
    return Cell(matrix, row, idx, mod);
}

Cell::Cell(const SparseMatrix *matrix_, size_t row, size_t column, bool mod_) throw() : tree(const_cast<SparseMatrix *>(matrix_)->tree), position(row, column), mod(mod_) {
    node = tree.get_node(position);
}

Cell& Cell::operator=(double rv) throw(std::logic_error) {
    if (mod == false) {
        throw std::logic_error("constant value could not be modified");
    }
    if (node == RBTree<MatrixIndex,double>::nil_node) {
        if (fabs(rv - 0.0) > SparseMatrix::get_precision()) {
            node = tree.insert(position, rv);
        }
    } else {
        if (fabs(rv - 0.0) > SparseMatrix::get_precision()) {
            node->item = rv;
        } else {
            tree.remove(node);
            node = RBTree<MatrixIndex,double>::nil_node;
        }
    }
    return *this;
}

Cell::operator double() const {
    if (node == RBTree<MatrixIndex,double>::nil_node) {
        return 0.0;
    }
    return node->item;
}

Cell& Cell::operator=(const Cell& rv) throw(std::logic_error) {
    if (this == &rv) {
        return *this;
    }
    return operator=(double(rv));
}

Cell& Cell::perform_operation(void (*op)(double&,double), double rv) throw(std::logic_error, std::bad_alloc) {
    bool node_exist = (node != RBTree<MatrixIndex,double>::nil_node);
    double tmp = (node_exist ? node->item : 0.0);

    if (mod == false) {
        throw std::logic_error("constant value could not be modified");
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
        node = RBTree<MatrixIndex,double>::nil_node;
    }
    return *this;
}

Cell& Cell::operator+=(double rv) throw(std::logic_error, std::bad_alloc) {
    return perform_operation([](double &lv, double rv_) {lv += rv_;}, rv);
}

Cell& Cell::operator-=(double rv) throw(std::logic_error, std::bad_alloc) {
    return perform_operation([](double &lv, double rv_) {lv -= rv_;}, rv);
}

Cell& Cell::operator*=(double rv) throw(std::logic_error) {
    return perform_operation([](double &lv, double rv_) {lv *= rv_;}, rv);
}

Cell& Cell::operator/=(double rv) throw(std::logic_error, std::domain_error) {
    if (fabs(rv - 0.0) <= SparseMatrix::get_precision()) {
        throw std::domain_error("could not divide by 0");
    }
    return perform_operation([](double &lv, double rv_) {lv /= rv_;}, rv);
}

bool Cell::operator==(double rv) const throw() {
    return (fabs(operator double() - rv) <= SparseMatrix::get_precision());
}

bool Cell::operator==(Cell& rv) const throw() {
    return (fabs(operator double() - double(rv)) <= SparseMatrix::get_precision());
}

bool Cell::operator!=(double rv) const throw() {
    return (fabs(operator double() - rv) > SparseMatrix::get_precision());
}

bool Cell::operator!=(Cell& rv) const throw() {
    return (fabs(operator double() - double(rv)) > SparseMatrix::get_precision());
}

template class Selector<ColumnSelector>;
template const Selector<ColumnSelector>& operator+(const Selector<ColumnSelector>&, size_t);
template const Selector<ColumnSelector>& operator+(size_t, const Selector<ColumnSelector>&);
template const Selector<ColumnSelector>& operator-(const Selector<ColumnSelector>&, size_t);

template class Selector<Cell>;
template const Selector<Cell>& operator+(const Selector<Cell>&, size_t);
template const Selector<Cell>& operator+(size_t, const Selector<Cell>&);
template const Selector<Cell>& operator-(const Selector<Cell>&, size_t);

