/*** Matrix.cc ***/

#include"matrix.h"
#include"myVector.h"
#include<iostream>
#include<iomanip>
#include<string>
#include<sstream>
#include<cmath>

std::ostream& operator<< (std::ostream& os, const Matrix& mat) {
	for(int i = 0; i < mat._row_size ; i++){
		os << "|";
		mat._row_list[i].printNobraket(os);
		os << "|"
			<< std::endl;
	}
	return os;
}

Matrix::Matrix(int n_row, int n_col) : _row_size(n_row), _col_size(n_col){
	// Row vector with all elements equal to 0
	RowVector temp(n_col);
	
	for(int i = 0; i < n_row; i++){
		_row_list.push_back(temp);
	}
}

Matrix::Matrix(std::vector<RowVector> row_vec) :_row_list(row_vec){
	int num_row = row_vec.size();

	if(num_row > 0){
		_row_size = num_row;
		_col_size = row_vec[0].size();
		for(int i = 1; i < num_row; i++){
			if(row_vec[i].size() != _col_size){
				std::cerr
					<< "ERROR (Matrix::Matrix): Number of columns differ among row vectors"
					<< std::endl;
				exit(1);
			}
		}
	} else{
		std::cerr 
			<< "ERROR (Matrix::Matrix): At least one row vector must be specified to construct a matrix"
			<< std::endl;
		exit(1);
	}
	
}

Matrix::Matrix(std::string str){
	if(str.empty() || str.at(0) != '{' || str.at(str.length() - 1) != '}'){
		std::cerr << "ERROR (Matrix::Matrix): Invalid input format"
			<< std::endl;
		exit(1);
	} 
	// Remove brackets {}
	str.erase(0,1);
	str.erase(str.length()-1,1);

	if(str.at(0) == ' '){str.erase(0,1);}
	if(str.at(str.length()-1) == ' '){str.erase(str.length()-1,1);}

	std::stringstream srm1(str);
	std::string row_str;
	std::vector<RowVector> row_vec;

	// Skip strings before the first { 
	std::getline(srm1, row_str, '{');
	while(std::getline(srm1, row_str, '}')){
		// Remove spaces at the beginning and the end
		if(row_str.at(0) == ' '){row_str.erase(0,1);}
		if(row_str.at(row_str.length()-1) == ' '){row_str.erase(str.length()-1,1);}

		row_vec.push_back(RowVector(row_str));

		// Skip strings before the next { 
		std::getline(srm1, row_str, '{');
	}

	(*this) = Matrix(row_vec);
}

void Matrix::setIdentity(){
	for(int i = 0; i < _row_size; i++) {
		_row_list[i].setOne(i);
	}
}

int Matrix::size() const{
	if(_row_size != _col_size){
		std::cerr << "ERROR (Matrix::size): The matrix is not square"
		<< std::endl;

		exit(1);
	}else{
		return _row_size;
	}
}

double Matrix::getMaxAbs() const{
	double max_abs = 0;	
	double temp_abs;

	for(int i = 0; i < _row_size; i++){
		temp_abs = _row_list[i].getMaxAbs();
		if(temp_abs > max_abs){max_abs = temp_abs;}
	}

	return max_abs;
}

double Matrix::getMaxAbsDiag() const{
	if(_row_size != _col_size){ 
		std::cerr << "ERROR (Matrix::getMaxDiag): The matrix is not square"
		<< std::endl;

		exit(1);
	}

	Matrix me = (*this);

	double max_abs_diag = 0;
	double temp_abs_diag;
	for(int i = 0; i < _row_size; i++){
		temp_abs_diag = fabs(me[i][i]);
		if( temp_abs_diag > max_abs_diag){
			max_abs_diag = temp_abs_diag;
		}
	}
	return max_abs_diag;
}

bool Matrix::isStochastic(double tolerance) const{
	ColVector e(_row_size);
	e.setAllOne();

	ColVector row_sum( (*this)*e );
	double min_sum = row_sum.getMinValue();
	double max_sum = row_sum.getMaxValue();

	if(max_sum <= 1 && min_sum > 1 - tolerance){return true;}
	else{return false;}
}

Matrix Matrix::getInverse() const{
	if(_row_size != _col_size){
		std::cerr << "ERROR (Matrix::getInverse): The matrix is not square"
		<< std::endl;

		exit(1);
	}
	int size = _row_size;

	Matrix me = (*this);

	Matrix inv(size, size);
	inv.setIdentity();

	std::vector<bool> fixed(size, false);

	while(1){
		int pivot = 0;
		double max_val = 0;

		for(int i=0; i < size; i++){
			if(fixed[i]){continue;}

			if(fabs(me[i][i]) > max_val){
				pivot = i;
				max_val = fabs(me[i][i]);
			}
		}

		if(max_val == 0){
			std::cerr << "ERROR (Matrix::getInverse): The matrix is not full rank"
			<< std::endl;

			exit(1);
		}

		double w = me[pivot][pivot];


		for(int j=0; j < size; j++){
			me[pivot][j] = me[pivot][j]/w;

			inv[pivot][j] = inv[pivot][j]/w;
		}

		for(int i=0; i < size; i++){
			if(i == pivot){continue;}


			double v = me[i][pivot];

			for(int j=0; j < size; j++){
				me[i][j] = me[i][j] - me[pivot][j]*v; 

				inv[i][j] = inv[i][j] - inv[pivot][j]*v;
			}
		}

		fixed[pivot] = true;

		bool finished = true;
		for(unsigned int i=0; i < fixed.size(); i++){
			finished = finished && fixed[i];
		}

		if(finished){break;}
	}
	
	return inv;

}

RowVector& Matrix::operator[](int k) {
	int num = _row_list.size();
	if(k < 0 || k >= num){
		std::cerr << "ERROR (Matrix::operator[]): The matrix does not have "<< k << "th row"
		<< " (It has [0--" << num - 1 << "]th row)"
		<< std::endl;

		exit(1);
	}

	return _row_list[k];
}

ColVector Matrix::operator*(const ColVector& col) const{
	std::vector<double> temp_vec;

	for(int i = 0; i < _row_size; i++){
		temp_vec.push_back(_row_list[i] * col);
	}
	return ColVector(temp_vec);
}

Matrix Matrix::operator+(const Matrix& mat) const{
	int num_row1= this->rowSize();
	int num_col1= this->colSize();
	int num_row2= mat.rowSize();
	int num_col2= mat.colSize();
	
	if(num_row1 != num_row2 || num_col1 != num_col2){
		std::cerr << "ERROR: Addition between inappropriate size matrices"
			<< std::endl;
			exit(1);
	}

	Matrix me = (*this);
	Matrix copy = mat;

	std::vector<RowVector> row_vec;
	for(int i = 0; i < num_row1; i++){
		row_vec.push_back(me[i] + copy[i]);
	}
	return Matrix(row_vec);
}

Matrix Matrix::operator*(const Matrix& mat) const{
	int num_row1= this->rowSize();
	int num_col1= this->colSize();
	int num_row2= mat.rowSize();

	if(num_col1 != num_row2){
		std::cerr << "ERROR: Multiplication between inappropriate size matrices"
			<< std::endl;
			exit(1);
	}

	std::vector<RowVector> row_vec;
	for(int i = 0; i < num_row1; i++){
		row_vec.push_back(this->_row_list[i] * mat);
	}
	return Matrix(row_vec);
}

Matrix Matrix::operator*(const double& c) const{
	std::vector<RowVector> new_vec;
	for(int i = 0; i < _row_size; i++){
		new_vec.push_back(_row_list[i] * c);
	}
	return Matrix(new_vec);
}

Matrix operator*(double c, const Matrix& mat){
	return mat*c;
}



void Matrix::makeStochastic(){
	ColVector e(_col_size);
	e.setAllOne();

	for(int i = 0; i < _row_size; i++){
		_row_list[i] = _row_list[i] * (1.0 / (_row_list[i] *e));
	}
}



