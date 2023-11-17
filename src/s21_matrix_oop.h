#ifndef SRC_S21_MATRIX_OOP_H
#define SRC_S21_MATRIX_OOP_H

#include "cmath"
#include "iostream"
#include "stdexcept"
#include "utility"

class S21Matrix {
 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&other) noexcept;
  ~S21Matrix();

  S21Matrix &operator=(const S21Matrix &other);
  S21Matrix &operator=(S21Matrix &&other) noexcept;

  int getRows() const;
  void setRows(int rows);
  int getCols() const;
  void setCols(int cols);

  S21Matrix &operator+=(const S21Matrix &other);
  S21Matrix &operator-=(const S21Matrix &other);
  S21Matrix &operator*=(const S21Matrix &other);
  S21Matrix &operator*=(const double num);

  S21Matrix operator+(const S21Matrix &other);
  S21Matrix operator-(const S21Matrix &other);
  S21Matrix operator*(const S21Matrix &other);
  S21Matrix operator*(const double num);

  bool operator==(const S21Matrix &other) const;
  double &operator()(int rows, int cols);

  bool EqMatrix(const S21Matrix &other) const;
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix &other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  void fillMatrix(const double *values);

 private:
  int rows_, cols_;
  double **matrix_;
  const double EPS = 1e-7;

  bool multValidation(const S21Matrix &other) const {
    return ((cols_ == other.rows_ && matrixOk(*this) && matrixOk(other)) ? 1
                                                                         : 0);
  }
  bool matrixOk(const S21Matrix &other) const {
    return ((other.rows_ > 0 && other.cols_ > 0 && other.matrix_ != nullptr)
                ? 1
                : 0);
  };
  bool sumSubValidation(const S21Matrix &other) const {
    return ((rows_ == other.rows_ && cols_ == other.cols_ && matrixOk(*this) &&
             matrixOk(other))
                ? 1
                : 0);
  }
  bool isSquareAndExist() const {
    return ((cols_ == rows_ && matrixOk(*this)) ? 1 : 0);
  }
  S21Matrix getMinor(const S21Matrix &other, int row, int column);
  void Free() noexcept;
  void memAlloc();
  void copyMatrix(const S21Matrix &other);
};

#endif  // SRC_S21_MATRIX_OOP_H
