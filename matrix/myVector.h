/*** MyVector.h ***/

#ifndef _MY_VECTOR_H
#define _MY_VECTOR_H

class BaseVector;
class RowVector;
class ColVector;

#include<vector>
#include<iostream>
#include<string>
#include<algorithm>
#include"matrix.h"


/* Base class for row and column vectors */
class BaseVector {
	protected:
		// List of elements
		std::vector<double> _element;
	public:
		// Constructor from the number of elements (set as zero)
		BaseVector(int num) : _element(num){this->setZero();}

		// Constructor from a list of elements
		BaseVector(std::vector<double> e){_element = e;}

		// Constructor from a string in the format {a_1, a_2, ..., a_n}
		BaseVector(std::string str);

		// Default constructor
		BaseVector(){}

		~BaseVector(){}

		// Sets all elements to be zero
		void setZero();
		// Sets the k-th element to be 1 and other elements to be 0
		void setOne(int k);

		// Returns the number of elements
		int size() const{return _element.size();}

		// Returns the set of elements as std::vector<double> 
		std::vector<double> getElements() const{return _element;}

		// Returns the minimum element
		double getMinValue() const;

		// Returns the maximum element
		double getMaxValue() const;

		// Returns the maximum absolute value
		double getMaxAbs() const{
			return std::max( fabs(this->getMinValue()), fabs(this->getMaxValue()) );
		}

		// Shift the elements to the right
		// (The last elements is removed and the first element becomes 0)
		void shift();

		double& operator[](int k);

		// Addition of vectors
		BaseVector operator+(const BaseVector& obj) const;

		// Subtraction of vectors
		BaseVector operator-(const BaseVector& obj) const;

		// Scalar multiplication
		BaseVector operator*(const double& c) const;

		// Scalar division
		BaseVector operator/(const double& c) const;

};

/* Class for row vectors */
class RowVector : public BaseVector {
	friend std::ostream& operator<<
			(std::ostream& os, const RowVector& row);
	public:
    // Constructor from the number of elements (set as zero)
		RowVector(int num) : BaseVector(num){}

		// Constructor from a list of elements
		RowVector(std::vector<double> e) : BaseVector(e){}

    // Constructor from a BaseVector
		RowVector(BaseVector b) : BaseVector(b){}

		// Default constructor
		RowVector() : BaseVector(){}
	
		// Print elements without brackets
		void printNobraket(std::ostream& os) const;

		// Print rounded elements without brackets
		void printIntNobraket(std::ostream& os) const;

		// Addition of vectors
			RowVector operator+(const RowVector& row) const{
				return RowVector(BaseVector::operator+(row));
			}
			void operator+(const BaseVector& b) const{
				std::cerr << "ERROR: Addition between vectors of different types" << std::endl;
				exit(1);
			}
		// Subtraction of vectors
			RowVector operator-(const RowVector& row) const{
				return RowVector(BaseVector::operator-(row));
			}
			void operator-(const BaseVector& b) const{
				std::cerr << "ERROR: Subtraction between vectors of different types" << std::endl;
				exit(1);
			}
		// Scalar multiplication
		RowVector operator*(const double& c) const{
			return RowVector(BaseVector::operator*(c));
		}

		// Scalar division
		RowVector operator/(const double& c) const{
			return RowVector(BaseVector::operator/(c));
		}

		// Inner product
			double operator*(const ColVector& col) const;
			void operator*(const BaseVector& b) const{
				std::cerr << "ERROR: Dot product between inappropriate type vectors" << std::endl;
				exit(1);
			}

		// Product with a matrix
		RowVector operator*(const Matrix& mat) const;

		// Comparison of vectors
		bool operator<(const RowVector& row) const;

		// Inner product for a given range
		double rangeProduct(int row_start, int row_end, const ColVector& col, int col_start, int col_end) const;
		double rangeProduct(int row_start, int row_end, const BaseVector& b, int col_start, int col_end) const{
				std::cerr << "ERROR: Range product between inappropriate type vectors" << std::endl;
				exit(1);
		}



};

/* Class for column vectors */
class ColVector : public BaseVector {
		// Necessary to define inner product 
		friend class RowVector;

		friend std::ostream& operator<<
			(std::ostream& os, const ColVector& col);

	public:
    // Constructor from the number of elements (set as zero)
		ColVector(int num) : BaseVector(num){}

    // Constructor from a list of elements
		ColVector(std::vector<double> e) : BaseVector(e){}

    // Constructor from a BaseVector 
		ColVector(BaseVector b) : BaseVector(b){}

		// Default constructor
		ColVector() : BaseVector(){}

		// Set all elements to be one
		void setAllOne(){
			int num = this->size();
			for (int i = 0; i < num; i++){_element[i] = 1;}
		}

		// Addition of vectors
			ColVector operator+(const ColVector& col) const{
				return ColVector(BaseVector::operator+(col));
			}
			void operator+(const BaseVector& b) const{
				std::cerr << "ERROR: Addition between vectors of different types" 
					<< std::endl;
				exit(1);
			}

		// Subtraction of vectors
			ColVector operator-(const ColVector& col) const{
				return ColVector(BaseVector::operator-(col));
			}
			void operator-(const BaseVector& b) const{
				std::cerr << "ERROR: Subtraction between vectors of different types" 
					<< std::endl;
				exit(1);
			}
	
		// Scalar multiplication
		ColVector operator*(const double& c) const{
			return ColVector(BaseVector::operator*(c));
		}

		// Scalar division
		ColVector operator/(const double& c) const{
			return ColVector(BaseVector::operator/(c));
		}

		// Multiply with a row vector from the left
			Matrix operator*(const RowVector& row) const;
			void operator*(const BaseVector& b) const{
				std::cerr << "ERROR: Multiplication between inappropriate type vectors"
					<< std::endl;
				exit(1);
			}
};

#endif
