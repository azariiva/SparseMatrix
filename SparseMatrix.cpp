#include "SparseMatrix.hpp"
#include <cmath>
#include <iostream>

void SparseMatrix::set_precision(double ueps) {
    SparseMatrix::eps = ueps;
}

double SparseMatrix::get_precision() {
    return SparseMatrix::eps;
}

SparseMatrix::SparseMatrix(size_t uheight, size_t uwidth) {
    row = 0;
    height = uheight;
    width = uwidth;
}

SparseMatrix::SparseMatrix(const SparseMatrix& src) : tree(src.tree) {
    row = 0;
    height = src.height;
    width = src.width;
}

SparseMatrix& SparseMatrix::operator=(const SparseMatrix& src) {
    height = src.height;
    width = src.width;
    tree = src.tree;
    return *this;
}

double SparseMatrix::get(size_t urow, size_t column) const {
    Node<MatrixIndex, double> *node;

    if (urow >= height || column >= width) {
        //TODO: raise error
        std::cerr << "Out of bounds" << std::endl;
        return 0.0;
    }
    node = tree.get_node(MatrixIndex(urow, column));
    if (node == RBTree<MatrixIndex, double>::nil_node) {
        return 0.0;
    }
    return node->item;
}

void SparseMatrix::set(size_t urow, size_t column, double val) {
    if (urow >= height || column >= width) {
        // TODO: raise error
        std::cerr << "Out of bounds" << std::endl;
        return ;
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

SparseArray SparseMatrix::operator[](size_t urow) {
     if (row + urow >= height) {
        // TODO: raise error
        std::cerr << "Out of bounds exception: height exceeded(" << row + urow << ")" << std::endl;
    }
    return SparseArray(this, row + urow);
}

SparseArray SparseMatrix::operator[](size_t urow) const {
    if (row + urow >= height) {
        // TODO: raise error
        std::cerr << "Out of bounds exception: height exceeded(" << row + urow << ")" << std::endl;
    }
    return SparseArray(const_cast<SparseMatrix *>(this), row + urow, false);
}

SparseArray SparseMatrix::operator*() {
    if (row >= height) {
        // TODO: raise error
        std::cerr << "Out of bounds exception: height exceeded(" << row << ")" << std::endl;
    }
    return SparseArray(this, row);
}

SparseArray SparseMatrix::operator*() const {
    if (row >= height) {
        // TODO: raise error
        std::cerr << "Out of bounds exception: height exceeded(" << row << ")" << std::endl;
    }
    return SparseArray(const_cast<SparseMatrix *>(this), row, false);
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
        // TODO: raise error
        std::cerr << "ValueError: operands could not be broadcast together" << std::endl;
    }
    for (size_t row_ = 0; row_ < height; row_++) {
        for (size_t column = 0; column < width; column++) {
            (operator[](row_)[column]).perform_operation(op, rv[row_][column]);
        }
    }
    return *this;
}

/*
    Все 4 нижестоящие функции - обёртки вокруг SparseArrayProxy::perform_operation()
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
        std::cerr << "ValueError: operands could not be broadcast together" << std::endl;
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

SparseArray::SparseArray(SparseMatrix *um, size_t urow, bool umodifyable) {
    m = um;
    m->row = 0;
    column = 0;
    row = urow;
    node = nullptr;
    modifyable = umodifyable;
}

SparseArray& SparseArray::operator+(size_t rv) {
    column += rv;
    return *this;
}

SparseArray& SparseArray::operator-(size_t rv) {
    column -= rv;
    return *this;
}

SparseArrayProxy SparseArray::operator[](size_t ucolumn) {
    column += ucolumn;
    if (column >= m->width) {
        // TODO: raise error
        std::cerr << "Out of bounds exception: width exceeded" << std::endl;
    }
    node = m->tree.get_node(MatrixIndex(row, column));
    return SparseArrayProxy(this);
}

SparseArrayProxy SparseArray::operator*() {
     if (column >= m->width) {
        // TODO: raise error
        std::cerr << "Out of bounds exception: width exceeded" << std::endl;
    }
    node = m->tree.get_node(MatrixIndex(row, column));
    return SparseArrayProxy(this);
}

SparseArrayProxy::SparseArrayProxy(SparseArray *ureal) {
    real = ureal;
}

SparseArrayProxy& SparseArrayProxy::operator=(double val) {
    if (real->modifyable == false) {
        // TODO: raise error
        std::cerr << "Can't modify constant matrix" << std::endl;
        return *this;
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

SparseArrayProxy& SparseArrayProxy::operator=(const SparseArrayProxy& rv) {
    return operator=(double(rv));
}

SparseArrayProxy::operator double() const {
    // std::cout << "(" << real->row << ' ' << real->column << ")";
    if (real->node == Node<MatrixIndex,double>::nil_node) {
        return 0.0;
    }
    return real->node->item;
}

/*
    Выполнить операцию вида <op>=rv
*/
SparseArrayProxy& SparseArrayProxy::perform_operation(void (*op)(double&,double), double rv) {
    bool node_exist = (real->node != Node<MatrixIndex,double>::nil_node);
    double tmp = (node_exist ? real->node->item : 0.0);

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
    Все нижестоящие 4 функции - обёртки вокруг SparseArrayProxy::perform_operation()
*/
SparseArrayProxy& SparseArrayProxy::operator+=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv += rv_;}, rv);
}

SparseArrayProxy& SparseArrayProxy::operator-=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv -= rv_;}, rv);
}

SparseArrayProxy& SparseArrayProxy::operator*=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv *= rv_;}, rv);
}

SparseArrayProxy& SparseArrayProxy::operator/=(double rv) {
    return perform_operation([](double &lv, double rv_) {lv /= rv_;}, rv);
}

bool SparseArrayProxy::operator==(double rv) {
    return (fabs(operator double() - rv) < SparseMatrix::eps);
}

bool SparseArrayProxy::operator==(SparseArrayProxy& rv) {
    return (fabs(operator double() - double(rv)) < SparseMatrix::eps);
}

bool SparseArrayProxy::operator!=(double rv) {
    return (fabs(operator double() - rv) >= SparseMatrix::eps);
}

bool SparseArrayProxy::operator!=(SparseArrayProxy& rv) {
    return (fabs(operator double() - double(rv)) >= SparseMatrix::eps);
}

double SparseMatrix::eps  = 1e-8;
