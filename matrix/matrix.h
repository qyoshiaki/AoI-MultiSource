/*** Matrix.h ***/

#ifndef _MATRIX_H
#define _MATRIX_H

class Matrix;

#include<vector>
#include<iostream>
#include<cmath>
#include"myVector.h"


class Matrix{
	friend std::ostream& operator<<
			(std::ostream& os, const Matrix& mat);

	// Necessary to define multiplication with row vectors
	friend class RowVector;

	private:
		// Elements of the matrix (represented as list of row vectors)
		std::vector<RowVector> _row_list;
		// Number of rows
		int _row_size;
		// Number of columuns
		int _col_size;
		
	public:
		// Constructor (initialized as a zero matrix)
		Matrix(int n_row, int n_col);

		// Constructor from a list of row vectors
		Matrix(std::vector<RowVector> row_vec);

		// Constructor from a string in the format {{a_11, ..., a_1n}, ..., {a_n1, ..., a_nn}}
		Matrix(std::string str);

		// Default constructor
		Matrix() : _row_size(0), _col_size(0){}

		~Matrix(){}
		
		// Sets the matrix to be an identity matrix
		void setIdentity();

		// Returns the number of rows
		int rowSize() const{return _row_size;}
	
		// Returns the number of columns
		int colSize() const{return _col_size;}

		// Returns the size of the matrix (aborted if it is not a square matrix)
		int size() const;

		// Returns the maximum absolute value among all elements
		double getMaxAbs() const;

		// Returns the maximum absolute value among diagonal elements
		double getMaxAbsDiag() const;

		// Returns if the matrix is stochastic (tolerance denotes the allowable error)
		bool isStochastic(double tolerance) const;

		// Returns the inverse matrix
		Matrix getInverse() const;

		RowVector& operator[](int k);

		// Product with column vectors
		ColVector operator*(const ColVector& col) const;

		// Addition of matrices
		Matrix operator+(const Matrix& mat) const;

		// Subtraction of matrices
		Matrix operator-(const Matrix& mat) const{
			return (*this) + (mat * (-1));
		}

		// Product of matrices
		Matrix operator*(const Matrix& mat) const;

		// Scalar multiplication
		Matrix operator*(const double& c) const;


		// Normalize each row to be stochastic
		void makeStochastic();

};


// Scalar multiplication
Matrix operator*(double c, const Matrix& mat);

#endif
