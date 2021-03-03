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

SparseMatrixProxy SparseMatrix::operator[](size_t row) {
    return SparseMatrixProxy(this, row, true);
}

SparseMatrixProxy SparseMatrix::operator[](size_t row) const {
    return SparseMatrixRow(const_cast<SparseMatrix *>(this), row, true, false);
}

SparseMatrixProxy SparseMatrix::operator*() {
    return SparseMatrixProxy(this, 0, true);
}

SparseMatrixProxy SparseMatrix::operator*() const {
    return SparseMatrixRow(const_cast<SparseMatrix *>(this), 0, false);
}

SparseMatrixProxy SparseMatrix::operator+(size_t row) {
    return SparseMatrixProxy(this, row);
}

SparseMatrixProxy SparseMatrix::operator+(size_t row) const {
    return SparseMatrixRow(const_cast<SparseMatrix *>(this), row, false, false);
}

SparseMatrixProxy SparseMatrix::operator-(size_t row) {
    return SparseMatrixProxy(this, -row);
}

SparseMatrixProxy SparseMatrix::operator-(size_t row) const {
    return SparseMatrixRow(const_cast<SparseMatrix *>(this), -row, false, false);
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
    Все 4 нижестоящие функции - обёртки вокруг SparseMatrix::perform_operation()
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

SparseMatrixProxy::SparseMatrixProxy(SparseMatrix *m_, size_t row_, bool deref_, bool modifyable_) {
    m = srcm;
    row = row_;
    column = 0;
    node = NULL;
    deref = deref_;
    modifyable = modifyable_;
    added = 0;
}

SparseMatrixProxy& SparseMatrixProxy::operator+(size_t rv) {
    if (deref == false) {
        row += rv;
    } else if (node == NULL) {
        column += rv;
    } else {
        added += rv;
    }
    return *this;
}

SparseMatrixProxy& SparseMatrixProxy::operator-(size_t rv) {
    if (deref == false) {
        row -= rv;
    } else if (node == NULL) {
        column -= rv;
    } else {
        added -= rv;
    }
    return *this;
}

SparseMatrixProxy& SparseMatrixProxy::operator[](size_t i) {
    if (deref == false) {
        row += i;
        if (row >= m->height) {
            throw std::out_of_range("SparseMatrixProxy::operator[] : out of bounds exception(height exceeded)");
        }
        deref = true;
    } else if (node == NULL) {
        column += i;
        if (column >= m->width) {
            throw std::out_of_range("SparseMatrixProxy::operator[] : out of bounds exception(width exceeded)");
        }
        node = m->tree.get_node(MatrixIndex(row, column));
    } else {
        throw std::logic_error("SparseMatrixProxy::operator[] : could not dereference value");
    }
    return *this;
}

SparseMatrixProxy& SparseMatrixProxy::operator*() {
    try {
        return this->operator[](0);
    } catch(std::logic_error&) {
        throw std::logic_error("SparseMatrixProxy::operator* : could not dereference value");
    }
}

SparseMatrixProxy& SparseMatrixProxy::operator=(double val) {
    if (node == NULL) {
        throw std::logic_error("SparseMatrixProxy::operator= : lvalue could not be assigned")
    }
    if (modifyable == false) {
        throw std::logic_error("SparseMatrixRowProxy::operator= : constant instance could not be modified");
    }
    if (node == Node<MatrixIndex,double>::nil_node) {
        if (fabs(val - 0.0) > SparseMatrix::eps) {
            node = m->tree.insert(MatrixIndex(row, column), val);
        }
    } else {
        if (fabs(val - 0.0) > SparseMatrix::eps) {
            node->item = val;
        } else {
            m->tree.remove(real->node);
            node = Node<MatrixIndex,double>::nil_node;
        }
    }
    return *this;
}

SparseMatrixProxy& SparseMatrixProxy::operator=(const SparseMatrixProxy& rv) {
    if (this == &rv) {
        return *this;
    }
    return operator=(double(rv));
}

SparseMatrixRowProxy::operator double() const {
    if (node == NULL) {
        throw std::logic_error("SparseMatrixRowProxy::operator double() : could not convert to double")
    }
    if (node == Node<MatrixIndex,double>::nil_node) {
        return 0.0 + (ssize_t)added;
    }
    return real->node->item + (ssize_t)added;
}

/*
    Выполнить операцию вида <op>=rv
*/
SparseMatrixProxy& SparseMatrixProxy::perform_operation(void (*op)(double&,double), double rv) {
    bool node_exist = (node != Node<MatrixIndex,double>::nil_node);
    double tmp = (node_exist ? node->item : 0.0);

    if (modifyable == false) {
        throw std::runtime_error("SparseMatrixProxy::perform_operation : constant instance could not be modified");
    }
    op(tmp, rv);
    if (fabs(tmp - 0.0) > SparseMatrix::eps) {
        if (node_exist) {
            node->item = tmp;
        } else {
            node = m->tree.insert(MatrixIndex(row, column), tmp);
        }
    } else if (node_exist) {
        m->tree.remove(node);
        node = Node<MatrixIndex,double>::nil_node;
    }
    return *this;
}

/*
    Все нижестоящие 4 функции - обёртки вокруг SparseMatrixRowProxy::perform_operation()
*/
SparseMatrixProxy& SparseMatrixProxy::operator+=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv += rv_;}, rv);
}

SparseMatrixProxy& SparseMatrixProxy::operator-=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv -= rv_;}, rv);
}

SparseMatrixProxy& SparseMatrixProxy::operator*=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv *= rv_;}, rv);
}

SparseMatrixProxy& SparseMatrixProxy::operator/=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv /= rv_;}, rv);
}

bool SparseMatrixProxy::operator==(double rv) const {
    return (fabs(operator double() - rv) < SparseMatrix::eps);
}

bool SparseMatrixProxy::operator==(SparseMatrixProxy& rv) const {
    return (fabs(operator double() - double(rv)) < SparseMatrix::eps);
}

bool SparseMatrixProxy::operator!=(double rv) const {
    return (fabs(operator double() - rv) >= SparseMatrix::eps);
}

bool SparseMatrixProxy::operator!=(SparseMatrixProxy& rv) const {
    return (fabs(operator double() - double(rv)) >= SparseMatrix::eps);
}

SparseMatrixProxy& operator+(size_t lv, const SparseMatrixProxy& rv) {
    return const_cast<SparseMatrixProxy&>(rv).operator+(lv);
}

SparseMatrix& operator+(size_t lv, SparseMatrix& rv) {
    return rv.operator+(lv);
}

const SparseMatrix& operator+(size_t lv, const SparseMatrix& rv) {
    return rv.operator+(lv);
}

double SparseMatrix::eps  = 1e-8;
size_t SparseMatrix::instance_quantity = 0;