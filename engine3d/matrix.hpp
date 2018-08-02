
#ifndef OCTILLION_ENGINE3D_MATRIX_HEADER
#define OCTILLION_ENGINE3D_MATRIX_HEADER

#include <valarray>

template <class T>
class Matrix
{
public:
	Matrix() {};

	Matrix(const Matrix& matrix)
	{
		row_ = matrix.row_;
		column_ = matrix.column_;
		matrix_ = matrix.matrix_;
	}

	Matrix& operator=(const Matrix& matrix)
	{
		row_ = matrix.row_;
		column_ = matrix.column_;
		matrix_ = matrix.matrix_;
		return *this;
	}

	// create a matrix by row and column
	Matrix(size_t row, size_t column) : row_(row), column_(column)
	{
		matrix_.resize(row_* column_);
	}

	// create a matrix with 1 row and 4 column (i.e. point matrix)
	Matrix(T x, T y, T z)
	{
		row_ = 1;
		column_ = 4;
		matrix_.resize(4);
		matrix_[0] = x;
		matrix_[1] = y;
		matrix_[2] = z;
		matrix_[3] = 1;
	}

	// create a matrix by row, column and a set of data
	Matrix(size_t row, size_t column, std::valarray<T> arr) : row_(row), column_(column)
	{
		matrix_ = arr;
	}

	// assign matrix value by row, column and an array of data
	void set(int row, int column, T data[])
	{
		row_ = row;
		column_ = column;
		set(data);			
	}

	// assign matrix value by an array of data without changing size
	void set(T data[])
	{
		std::valarray<T> matrix(data, row_*column_);
		matrix_ = matrix;
	}

	// member data accessor, row size
	size_t row() const { return row_; }

	// member data accessor, column size
	size_t column() const { return column_; }

	// member data accessor, get value
	T at(size_t row, size_t column) const
	{
		return matrix_[row*row_ + column];
	}

	// member data accessor, get value
	T at(size_t index) const
	{
		return matrix_[index];
	}

	// member data accessor, set value
	void set(size_t row, size_t column, T val)
	{
		matrix_[row*row_ + column] = val;
	}

	// operator overloading, matrixA = matrixB * matrixC
	template <class T> 
	friend Matrix<T> operator* (const Matrix<T>& lhs, const Matrix<T>& rhs)
	{
		std::valarray<T> result(lhs.row_ * lhs.column_);

		for (size_t row = 0; row < lhs.row_; row++)
		{
			for (size_t column = 0; column < lhs.column_; column++)
			{
				std::valarray<T> lhval = (lhs.matrix_)[std::slice(row*lhs.column_, lhs.column_, 1)];
				std::valarray<T> rhval = (rhs.matrix_)[std::slice(column, rhs.row_, rhs.column_)];
				std::valarray<T> multi = lhval * rhval;

				result[row*lhs.row_ + column] = multi.sum();
			}
		}

		Matrix<T> ret(lhs.row_, lhs.column_, result);
		return ret;
	}

	template <class T>
	friend Matrix<T> operator- (const Matrix<T>& lhs, const Matrix<T>& rhs)
	{
		std::valarray<T> result = lhs.matrix_ - rhs.matrix_;
		Matrix<T> ret(lhs.row_, lhs.column_, result);
		return ret;
	}

	// operator overloading, matrixA *= matrixB
	Matrix<T>& operator*= (const Matrix<T>& rhs)
	{
		std::valarray<T> result(row_ * column_);

		for (size_t row = 0; row < row_; row++)
		{
			for (size_t column = 0; column < column_; column++)
			{
				std::valarray<T> lhval = (matrix_)[std::slice(row*column_, column_, 1)];
				std::valarray<T> rhval = (rhs.matrix_)[std::slice(column, rhs.row_, rhs.column_)];
				std::valarray<T> multi = lhval * rhval;

				result[row*row_ + column] = multi.sum();
			}
		}

		matrix_ = result;

		return *this;
	}

protected:
	std::valarray<T> matrix_;
	size_t row_;
	size_t column_;
};

#endif