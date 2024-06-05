#include "matrix.h"

s21::Matrix::Matrix() : rows_(1), cols_(1) { CreateMatrix(); }

s21::Matrix::Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  CreateMatrix();
}

s21::Matrix::Matrix(const s21::Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  CreateMatrix();
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

s21::Matrix::Matrix(s21::Matrix&& other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

s21::Matrix::~Matrix() { RemoveMatrix(); }

int s21::Matrix::GetRows() const { return rows_; }

int s21::Matrix::GetCols() const { return cols_; }

void s21::Matrix::SetRows(int rows) {
  if (rows <= 0) {
    throw std::out_of_range("Incorrect input. Rows must be greater than 0");
  }
  if (rows != rows_) {
    s21::Matrix result(rows, cols_);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols_; j++) {
        if (i < rows_) {
          result.matrix_[i][j] = matrix_[i][j];
        }
      }
    }
    RemoveMatrix();
    *this = std::move(result);
  }
}

void s21::Matrix::SetCols(int cols) {
  if (cols <= 0) {
    throw std::out_of_range("Incorrect input. Cols must be greater than 0");
  }
  if (cols != cols_) {
    s21::Matrix result(rows_, cols);
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols; j++) {
        if (j < cols_) {
          result.matrix_[i][j] = matrix_[i][j];
        }
      }
    }
    RemoveMatrix();
    *this = std::move(result);
  }
}
double s21::Matrix::GetMinValue() {
  double min = matrix_[0][0];
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++)
      if (matrix_[i][j] < min) min = matrix_[i][j];
  return min;
}

double s21::Matrix::GetMaxValue() {
  double max = matrix_[0][0];
  for (int i = 0; i < rows_; i++)
    for (int j = 0; j < cols_; j++)
      if (matrix_[i][j] > max) max = matrix_[i][j];
  return max;
}

void s21::Matrix::SetValue(double value) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] = value;
    }
  }
}

void s21::Matrix::SetValue(std::deque<double>& values) {
  if (rows_ * cols_ != static_cast<int>(values.size())) {
    throw std::logic_error("Different size of matrix ande deque");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) matrix_[i][j] = values[i * cols_ + j];
  }
}

void s21::Matrix::SetRandomValue(int max_value) {
  srand(time(NULL));
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (max_value == 0)
        matrix_[i][j] = arc4random() % (101) / 100.0;
      else
        matrix_[i][j] = arc4random() % (max_value + 1);
    }
  }
}

void s21::Matrix::Print() {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      std::cout << matrix_[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

bool s21::Matrix::EqMatrix(const s21::Matrix& other) const noexcept {
  bool res = false;
  if (rows_ > 0 && cols_ > 0 && rows_ == other.rows_ && cols_ == other.cols_ &&
      matrix_ != nullptr && other.matrix_ != nullptr) {
    res = true;
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++)
        if (fabs(matrix_[i][j] - other.matrix_[i][j]) >= 1e-7) {
          res = false;
        }
      if (res == false) break;
    }
  }
  return res;
}

void s21::Matrix::SumMatrix(const s21::Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error("Different size of matrix");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void s21::Matrix::SubMatrix(const s21::Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error("Different size of matrix");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void s21::Matrix::MulNumber(const double num) noexcept {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] *= num;
    }
  }
}
void s21::Matrix::SumNumber(const double num) noexcept {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] += num;
    }
  }
}

void s21::Matrix::MulMatrix(const s21::Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::logic_error("Incorrect matrix sizes");
  }
  s21::Matrix result(rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      for (int k = 0; k < cols_; k++) {
        result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  RemoveMatrix();
  *this = std::move(result);
}

void s21::Matrix::Hadamard(const s21::Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::logic_error("Different size of matrix");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] *= other.matrix_[i][j];
    }
  }
}

s21::Matrix s21::Matrix::Transpose() noexcept {
  s21::Matrix result(cols_, rows_);
  for (int i = 0; i < cols_; ++i) {
    for (int j = 0; j < rows_; ++j) {
      result.matrix_[i][j] = matrix_[j][i];
    }
  }
  return result;
}

s21::Matrix s21::Matrix::operator+(const s21::Matrix& other) const {
  s21::Matrix result(*this);
  result.SumMatrix(other);
  return result;
}

s21::Matrix s21::Matrix::operator-(const s21::Matrix& other) const {
  s21::Matrix result(*this);
  result.SubMatrix(other);
  return result;
}

s21::Matrix s21::Matrix::operator*(const s21::Matrix& other) const {
  s21::Matrix result(*this);
  result.MulMatrix(other);
  return result;
}

s21::Matrix s21::Matrix::operator*(const double multiplicator) const noexcept {
  s21::Matrix result(*this);
  result.MulNumber(multiplicator);
  return result;
}

s21::Matrix operator*(const double multiplicator, const s21::Matrix& other) {
  s21::Matrix result = other;
  result.MulNumber(multiplicator);
  return result;
}

bool s21::Matrix::operator==(const s21::Matrix& other) const noexcept {
  return EqMatrix(other);
}

s21::Matrix& s21::Matrix::operator=(const s21::Matrix& other) noexcept {
  if (this != &other) {
    RemoveMatrix();
    rows_ = other.rows_;
    cols_ = other.cols_;
    CreateMatrix();
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
  return *this;
}

s21::Matrix& s21::Matrix::operator=(s21::Matrix&& other) noexcept {
  if (this != &other) {
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    std::swap(matrix_, other.matrix_);
  }
  return *this;
}

s21::Matrix& s21::Matrix::operator+=(const s21::Matrix& other) {
  SumMatrix(other);
  return *this;
}

s21::Matrix& s21::Matrix::operator+=(const double num) noexcept {
  SumNumber(num);
  return *this;
}

s21::Matrix& s21::Matrix::operator-=(const s21::Matrix& other) {
  SubMatrix(other);
  return *this;
}

s21::Matrix& s21::Matrix::operator*=(const s21::Matrix& other) {
  MulMatrix(other);
  return *this;
}

s21::Matrix& s21::Matrix::operator*=(const double multiplicator) noexcept {
  MulNumber(multiplicator);
  return *this;
}

double& s21::Matrix::operator()(int i, int j) {
  if ((rows_ <= i) || (cols_ <= j) || (i < 0) || (j < 0)) {
    throw std::out_of_range("Index out of range");
  }
  return matrix_[i][j];
}

s21::Matrix& s21::Matrix::operator%=(const s21::Matrix& other) {
  Hadamard(other);
  return *this;
}

double& s21::Matrix::operator()(const int i, const int j) const {
  if ((rows_ <= i) || (cols_ <= j) || (i < 0) || (j < 0)) {
    throw std::out_of_range("Index out of range");
  }
  return matrix_[i][j];
}

void s21::Matrix::CreateMatrix() {
  if ((rows_ <= 0) || (cols_ <= 0)) {
    throw std::out_of_range(
        "Incorrect input. Rows and colomns must be greater than 0");
  }
  matrix_ = new double*[rows_]();
  if (!matrix_) throw std::bad_alloc();
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
    if (!matrix_[i]) throw std::bad_alloc();
  }
}

void s21::Matrix::RemoveMatrix() {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
    matrix_ = nullptr;
  }
  rows_ = 0;
  cols_ = 0;
}
