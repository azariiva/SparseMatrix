#include "SparseMatrix.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

void SparseMatrix::set_precision(double ueps) {
    if (SparseMatrix::instance_quantity != 0) {
        throw std::runtime_error("SparseMatrix::set_precision : precision could not be modified when instances of SparseMatrix exist");
    } else {
        SparseMatrix::eps = ueps;
    }
}

double SparseMatrix::get_precision() {
    return SparseMatrix::eps;
}

SparseMatrix::SparseMatrix(size_t uheight, size_t uwidth) {
    row = 0;
    height = uheight;
    width = uwidth;
    SparseMatrix::instance_quantity++;
}

SparseMatrix::SparseMatrix(const SparseMatrix& src) : tree(src.tree) {
    row = 0;
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

double SparseMatrix::get(size_t urow, size_t column) const {
    Node<MatrixIndex, double> *node;

    if (urow >= height || column >= width) {
        throw std::out_of_range("SparseMatrix::get : out of bounds exception(height or width exceeded)");
    }
    node = tree.get_node(MatrixIndex(urow, column));
    if (node == RBTree<MatrixIndex, double>::nil_node) {
        return 0.0;
    }
    return node->item;
}

void SparseMatrix::set(size_t urow, size_t column, double val) {
    if (urow >= height || column >= width) {
        throw std::out_of_range("SparseMatrix::set : out of bounds exception(height or width exceeded)");
    }
    if (fabs(val - 0.0) <= SparseMatrix::eps) {
        tree.remove(MatrixIndex(urow, column));
    } else {
        tree.insert(MatrixIndex(urow, column), val);
    }
}

size_t SparseMatrix::num_rows() const {
    return height;
}

size_t SparseMatrix::num_columns() const {
    return width;
}

SparseMatrixRow SparseMatrix::operator[](size_t urow) {
     if (row + urow >= height) {
        throw std::out_of_range("SparseMatrix::operator[] : out of bounds exception(height exceeded)");
    }
    return SparseMatrixRow(this, row + urow);
}

SparseMatrixRow SparseMatrix::operator[](size_t urow) const {
    if (row + urow >= height) {
        throw std::out_of_range("SparseMatrix::operator[] : out of bounds exception(height exceeded)");
    }
    return SparseMatrixRow(const_cast<SparseMatrix *>(this), row + urow, false);
}

SparseMatrixRow SparseMatrix::operator*() {
    if (row >= height) {
        throw std::out_of_range("SparseMatrix::operator* : out of bounds exception(height exceeded)");
    }
    return SparseMatrixRow(this, row);
}

SparseMatrixRow SparseMatrix::operator*() const {
    if (row >= height) {
        throw std::out_of_range("SparseMatrix::operator* : out of bounds exception(height exceeded)");
    }
    return SparseMatrixRow(const_cast<SparseMatrix *>(this), row, false);
}

SparseMatrix& SparseMatrix::operator+(size_t urow) {
    row += urow;
    return *this;
}

const SparseMatrix& SparseMatrix::operator+(size_t urow) const {
    const_cast<SparseMatrix *>(this)->row += urow;
    return *this;
}

SparseMatrix& SparseMatrix::operator-(size_t urow) {
    row -= urow;
    return *this;
}

const SparseMatrix& SparseMatrix::operator-(size_t urow) const {
    const_cast<SparseMatrix *>(this)->row -= urow;
    return *this;
}

SparseMatrix& SparseMatrix::perform_operation(void (*op)(double&,double), const SparseMatrix&rv) {
    // TODO: optimize a little
    if (height != rv.height || width != rv.width) {
        throw std::invalid_argument("SparseMatrix::perform_operation : operands could not be broadcast together");
    }
    for (size_t row_ = 0; row_ < height; row_++) {
        for (size_t column = 0; column < width; column++) {
            (operator[](row_)[column]).perform_operation(op, rv[row_][column]);
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

SparseMatrixRow::SparseMatrixRow(SparseMatrix *um, size_t urow, bool umodifyable) {
    m = um;
    m->row = 0;
    column = 0;
    row = urow;
    node = nullptr;
    modifyable = umodifyable;
}

SparseMatrixRow& SparseMatrixRow::operator+(size_t rv) {
    column += rv;
    return *this;
}

SparseMatrixRow& SparseMatrixRow::operator-(size_t rv) {
    column -= rv;
    return *this;
}

SparseMatrixRowProxy SparseMatrixRow::operator[](size_t ucolumn) {
    column += ucolumn;
    if (column >= m->width) {
        throw std::out_of_range("SparseMatrixRow::operator[] : out of bounds exception(width exceeded)");
    }
    node = m->tree.get_node(MatrixIndex(row, column));
    return SparseMatrixRowProxy(this);
}

SparseMatrixRowProxy SparseMatrixRow::operator*() {
     if (column >= m->width) {
        throw std::out_of_range("SparseMatrixRow::operator* : out of bounds exception(width exceeded)");
    }
    node = m->tree.get_node(MatrixIndex(row, column));
    return SparseMatrixRowProxy(this);
}

SparseMatrixRowProxy::SparseMatrixRowProxy(SparseMatrixRow *ureal) {
    real = ureal;
}

SparseMatrixRowProxy& SparseMatrixRowProxy::operator=(double val) {
    if (real->modifyable == false) {
        throw std::runtime_error("SparseMatrixRowProxy::operator= : constant instance could not be modified");
    }
    if (real->node == Node<MatrixIndex,double>::nil_node) {
        if (fabs(val - 0.0) > SparseMatrix::eps) {
            real->node = real->m->tree.insert(MatrixIndex(real->row, real->column), val);
        }
    } else {
        if (fabs(val - 0.0) > SparseMatrix::eps) {
            real->node->item = val;
        } else {
            real->m->tree.remove(real->node);
            real->node = Node<MatrixIndex,double>::nil_node;
        }
    }
    return *this;
}

SparseMatrixRowProxy& SparseMatrixRowProxy::operator=(const SparseMatrixRowProxy& rv) {
    if (this == &rv) {
        return *this;
    }
    return operator=(double(rv));
}

SparseMatrixRowProxy::operator double() const {
    if (real->node == Node<MatrixIndex,double>::nil_node) {
        return 0.0;
    }
    return real->node->item;
}

/*
    Выполнить операцию вида <op>=rv
*/
SparseMatrixRowProxy& SparseMatrixRowProxy::perform_operation(void (*op)(double&,double), double rv) {
    bool node_exist = (real->node != Node<MatrixIndex,double>::nil_node);
    double tmp = (node_exist ? real->node->item : 0.0);

    if (real->modifyable == false) {
        throw std::runtime_error("SparseMatrixRowProxy::perform_operation : constant instance could not be modified");
    }
    op(tmp, rv);
    if (fabs(tmp - 0.0) > SparseMatrix::eps) {
        if (node_exist) {
            real->node->item = tmp;
        } else {
            real->node = real->m->tree.insert(MatrixIndex(real->row, real->column), tmp);
        }
    } else if (node_exist) {
        real->m->tree.remove(real->node);
        real->node = Node<MatrixIndex,double>::nil_node;
    }
    return *this;
}

/*
    Все нижестоящие 4 функции - обёртки вокруг SparseMatrixRowProxy::perform_operation()
*/
SparseMatrixRowProxy& SparseMatrixRowProxy::operator+=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv += rv_;}, rv);
}

SparseMatrixRowProxy& SparseMatrixRowProxy::operator-=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv -= rv_;}, rv);
}

SparseMatrixRowProxy& SparseMatrixRowProxy::operator*=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv *= rv_;}, rv);
}

SparseMatrixRowProxy& SparseMatrixRowProxy::operator/=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv /= rv_;}, rv);
}

bool SparseMatrixRowProxy::operator==(double rv) const {
    return (fabs(operator double() - rv) < SparseMatrix::eps);
}

bool SparseMatrixRowProxy::operator==(SparseMatrixRowProxy& rv) const {
    return (fabs(operator double() - double(rv)) < SparseMatrix::eps);
}

bool SparseMatrixRowProxy::operator!=(double rv) const {
    return (fabs(operator double() - rv) >= SparseMatrix::eps);
}

bool SparseMatrixRowProxy::operator!=(SparseMatrixRowProxy& rv) const {
    return (fabs(operator double() - double(rv)) >= SparseMatrix::eps);
}

SparseMatrixRow& operator+(size_t lv, const SparseMatrixRow& rv) {
    return const_cast<SparseMatrixRow&>(rv).operator+(lv);
}

SparseMatrix& operator+(size_t lv, SparseMatrix& rv) {
    return rv.operator+(lv);
}

const SparseMatrix& operator+(size_t lv, const SparseMatrix& rv) {
    return rv.operator+(lv);
}

SparseMatrix operator*(double lv, const SparseMatrix& rv) {
    return rv.operator*(lv);
}

double SparseMatrix::eps  = 1e-8;
size_t SparseMatrix::instance_quantity = 0;