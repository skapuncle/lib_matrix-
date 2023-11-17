#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(0) {}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  memAlloc();
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(nullptr) {
  memAlloc();
  copyMatrix(other);
}

S21Matrix::S21Matrix(S21Matrix &&other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() { Free(); }

// Matrix func

bool S21Matrix::EqMatrix(const S21Matrix &other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return false;
  }
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++)
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) > EPS) {
        return false;
      }

  return true;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (!sumSubValidation(other)) {
    throw std::logic_error("Matrixes have different size or uninitialized");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (!sumSubValidation(other)) {
    throw std::logic_error("Matrixes have different size or uninitialized");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  if (!matrixOk(*this)) {
    throw std::invalid_argument("Invalid matrix size");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (!multValidation(other)) {
    throw std::logic_error(
        "Rows of 1st matrix are not equal columns of the 2nd one");
  }
  S21Matrix temp(rows_, other.cols_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < other.cols_; ++j) {
      for (int k = 0; k < cols_; ++k) {
        temp.matrix_[i][j] += (matrix_[i][k] * other.matrix_[k][j]);
      }
    }
  }
  *this = std::move(temp);
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix result(cols_, rows_);
  for (int i = 0; i < result.rows_; ++i) {
    for (int j = 0; j < result.cols_; ++j) {
      result.matrix_[i][j] = matrix_[j][i];
    }
  }
  return result;
}

S21Matrix S21Matrix::CalcComplements() {
  if (!isSquareAndExist()) {
    throw std::logic_error("Not square matrix, cannot calculate");
  }
  S21Matrix result(rows_, cols_);
  if (rows_ == 1) {
    result.matrix_[0][0] = matrix_[0][0];
  } else {
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        S21Matrix temp = getMinor(*this, i, j);
        result.matrix_[i][j] = pow(-1, i + j) * temp.Determinant();
      }
    }
  }
  return result;
}

double S21Matrix::Determinant() {
  if (!isSquareAndExist()) {
    throw std::logic_error("Not square matrix, cannot calculate");
  }
  double det = 0;
  if (rows_ == 1) {
    det = matrix_[0][0];
  } else if (rows_ == 2) {
    det = matrix_[0][0] * matrix_[1][1] - matrix_[1][0] * matrix_[0][1];
  } else {
    for (int i = 0; i < cols_; ++i) {
      S21Matrix temp = getMinor(*this, i, 0);
      det += pow(-1, i + 2) * matrix_[i][0] * temp.Determinant();
    }
  }
  return det;
}

S21Matrix S21Matrix::getMinor(const S21Matrix &other, int row, int column) {
  S21Matrix result(other.rows_ - 1, other.cols_ - 1);
  for (int i = 0, ii = 0; i < other.rows_; ++i) {
    if (i == row) {
      continue;
    }
    for (int j = 0, jj = 0; j < other.cols_; ++j) {
      if (j == column) {
        continue;
      }
      result.matrix_[ii][jj] = other.matrix_[i][j];
      jj++;
    }
    ii++;
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  if (!matrixOk(*this)) {
    throw std::invalid_argument("Invalid matrix size");
  }
  if (!isSquareAndExist()) {
    throw std::logic_error("Not square matrix, cannot calculate");
  }
  double det = Determinant();
  S21Matrix result(rows_, cols_);
  if (fabs(det) < EPS) {
    throw std::logic_error("Determinant is 0");
  }
  if (rows_ == 1) {
    result.matrix_[0][0] = 1 / det;
  } else {
    result = CalcComplements();
    result = result.Transpose();
    result.MulNumber(1.0 / det);
  }
  return result;
}

// getters and setters

int S21Matrix::getRows() const { return rows_; }

int S21Matrix::getCols() const { return cols_; }

void S21Matrix::setRows(int rows) {
  if (rows > 0 && rows != rows_) {
    S21Matrix temp(rows, cols_);
    temp.copyMatrix(*this);
    *this = temp;
  }
}

void S21Matrix::setCols(int cols) {
  if (cols > 0 && cols != cols_) {
    S21Matrix temp(rows_, cols);
    temp.copyMatrix(*this);
    *this = temp;
  }
}

// operators

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  S21Matrix copy{other};
  *this = std::move(copy);
  return *this;
}

S21Matrix &S21Matrix::operator=(S21Matrix &&other) noexcept {
  if (this != &other) {
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    std::swap(matrix_, other.matrix_);
  }

  return *this;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  MulMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix result(*this);
  result.MulNumber(num);
  return result;
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  SumMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  SubMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

bool S21Matrix::operator==(const S21Matrix &other) const {
  return EqMatrix(other);
}

double &S21Matrix::operator()(int rows, int cols) {
  if (rows > rows_ || cols > cols_ || rows < 0 || cols < 0) {
    throw std::out_of_range("Out of range");
  }
  return matrix_[rows][cols];
}

// helpers

void S21Matrix::fillMatrix(const double *values) {
  for (int i = 0, el = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++, el++) {
      matrix_[i][j] = values[el];
    }
  }
}

void S21Matrix::memAlloc() {
  if (rows_ < 1 || cols_ < 1) {
    throw std::invalid_argument("Invalid matrix size");
  }
  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; ++i) {
    matrix_[i] = new double[cols_]();
  }
}

void S21Matrix::copyMatrix(const S21Matrix &other) {
  int r_end = rows_ <= other.rows_ ? rows_ : other.rows_;
  size_t shift = cols_ <= other.cols_ ? cols_ : other.cols_;
  for (int i = 0; i < r_end; ++i) {
    std::copy(other.matrix_[i], other.matrix_[i] + shift, matrix_[i]);
  }
}

void S21Matrix::Free() noexcept {
  delete[] matrix_;
  rows_ = 0;
  cols_ = 0;

  matrix_ = nullptr;
}