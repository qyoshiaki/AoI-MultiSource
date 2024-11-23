/*** MyVector.cc ***/

#include"myVector.h"
#include<cstdlib>
#include<iostream>
#include<iomanip>
#include<string>
#include<sstream>
#include<vector>

BaseVector::BaseVector(std::string str){
	if(str.empty()){
		std::cerr << "ERROR (BaseVector::BaseVector): Invalid input format"
			<< std::endl;
		exit(1);
	} 
	// Remove {}
	if(str.at(0) == '{'){str.erase(0,1);}
	if(str.at(str.length()-1) == '}'){str.erase(str.length()-1,1);}

	std::stringstream srm1(str);
	std::string value_str;

	while(std::getline(srm1, value_str, ',')){
		std::stringstream srm2(value_str);

		double value;

		srm2 >> value;
		_element.push_back(value);
	}
}


void BaseVector::setZero(){
	int num = _element.size();
	for(int i = 0; i < num; i++){
		_element[i] = 0;
	}
}

void BaseVector::setOne(int k){
	int num = this->size();

	if (k < 0 || k  >= num){
		std::cerr 
			<< "ERROR (BaseVector::setOne): The vector does not have "<< k << "th element"
			<< " (It has [0--" << num - 1 << "]th element)"
			<< std::endl;
		exit(1);
	}
	else {
		for (int i = 0; i < num; i++){
			_element[i] = 0;
		}
		_element[k] = 1;
	}
}


double BaseVector::getMinValue() const{
	double min_value = _element[0];

	int num = _element.size();
	for(int i = 1; i < num; i++){
		if(_element[i] < min_value){
			min_value = _element[i];
		}
	}
	return min_value;
}

double BaseVector::getMaxValue() const{
	double max_value = _element[0];

	int num = _element.size();
	for(int i = 1; i < num; i++){
		if(_element[i] > max_value){
			max_value = _element[i];
		}
	}
	return max_value;
}

void BaseVector::shift(){
	if(_element.size() == 0){return;}

	for(int i=_element.size()-1; i > 0; i--){
		_element[i] = _element[i-1];
	}
	_element[0] = 0;
}

double& BaseVector::operator[](int k){

	int num = _element.size();
	if(k < 0 || k >= num){
		std::cerr << "ERROR (RowVector::operator[]): The vector does not have "<< k << "th element"
			<< " (It has [0--" << num - 1 << "]th element)"
			<< std::endl;
	}

	return _element[k];
}


BaseVector BaseVector::operator+(const BaseVector& obj) const{
	int num1 = this->size();
	int num2 = obj.size();

	if(num1 == num2){

		std::vector<double> temp(num1);

		for(int i = 0; i < num1; i++){
			temp[i] = this->_element[i] + obj._element[i];
		}
		return BaseVector(temp);
	} else{
		std::cerr 
			<< "ERROR (BaseVector::operator+): Addition between different size vectors"
			<< std::endl;
		exit(1);
	}
}

BaseVector BaseVector::operator-(const BaseVector& obj) const{
	int num1 = this->size();
	int num2 = obj.size();

	if(num1 == num2){

		std::vector<double> temp(num1);

		for(int i = 0; i < num1; i++){
			temp[i] = this->_element[i] - obj._element[i];
		}
		return BaseVector(temp);
	} else{
		std::cerr 
			<< "ERROR (BaseVector::operator-): Subtraction between different size vectors"
			<< std::endl;
		exit(1);
	}
}

BaseVector BaseVector::operator*(const double& c) const{
	int num = this->size();
	std::vector<double> value;

	for(int i = 0; i < num; i++){
		value.push_back(_element[i] * c);
	}
	return BaseVector(value);
}

BaseVector BaseVector::operator/(const double& c) const{
  if(c == 0){
    std::cerr << "ERROR (BaseVector::operator/): Zero division" << std::endl;
  }
  return (*this)*(1.0/c);
}

std::ostream& operator<<(std::ostream& os, const RowVector& row){
	int num = row.size();

	os << "[" ;

	for(int i = 0; i < num - 1; i++){
		os << row._element[i] << ", ";
	}

	os << row._element[num - 1] << "]";
	return os;
}

void RowVector::printNobraket(std::ostream& os) const{
	int num = this->size();

	for(int i = 0; i < num-1; i++){
		os << _element[i] << "  ";
	}
	os << _element[num-1];
}

void RowVector::printIntNobraket(std::ostream& os) const{
	int num = this->size();

	for(int i = 0; i < num-1; i++){
		os << (int)(_element[i]) << ", ";
	}
	os << (int)(_element[num-1]);

}


double RowVector::operator*(const ColVector& col) const{
	int num1 = this->size();
	int num2 = col.size();

	if(num1 == num2){
		double dProduct = 0;
		
		for(int i = 0; i < num1; i++){
			dProduct += this->_element[i] * col._element[i];
		}
		return dProduct;
	} else{
		std::cerr 
			<< "ERROR (RowVector::operator*): Dot product between different size vectors"
			<< std::endl;
		exit(1);
	}
}

double RowVector::rangeProduct(int row_start, int row_end, const ColVector& col, int col_start, int col_end) const{
	
	int num1 = row_end - row_start + 1;
	int num2 = col_end - col_start + 1;

	if(row_start < 0 || row_end >= (int)_element.size() || col_start < 0 || col_end >= (int)col._element.size() || num1 != num2){
		std::cerr 
			<< "ERROR (RowVector::rangeProduct): Entered range is Invalid "
			<< std::endl;
		exit(1);
	}

	double dProduct = 0;

	for(int i=0; i < num1; i++){
		dProduct += _element[row_start+i] * col._element[col_start+i];
	}

	return dProduct;

}




RowVector RowVector::operator*(const Matrix& mat) const{
	int num = this->size();

	if(mat.rowSize() != num){
		std::cerr
			<< "ERROR (RowVector::operator*): Multiplication between vector and matrix with different sizes"
			<< std::endl;
		exit(1);
	}

	RowVector new_row(mat.colSize());
	
	for(int i = 0; i < num; i++){
		new_row = new_row + mat._row_list[i] *  _element[i];
	}
	return new_row;
}

bool RowVector::operator<(const RowVector& row) const{
	if(this->size() != row.size()){
		std::cerr
			<< "ERROR (RowVector::operator<): Comparison between different size vectors"
			<< std::endl;
		exit(1);
	}

	for(int i = this->size() - 1; i >= 0; i--){
		if(_element[i] != row._element[i]){
			return (_element[i] < row._element[i]);
		}
	}
	return false;
}
	

std::ostream& operator<< (std::ostream& os, const ColVector& col){
	int num = col.size();

	for(int i = 0; i < num ; i++){
		os << "|" 
			<< std::setw(4) 
			<< col._element[i] 
			<< "|"
			<< std::endl;
	}
	return os;
}

Matrix ColVector::operator*(const RowVector& row) const{
	int num_row = this->size();
	std::vector<RowVector> row_vec;
		
	for(int i = 0; i < num_row; i++){
		row_vec.push_back(row * this->_element[i]);
	}
	return Matrix(row_vec);
}
