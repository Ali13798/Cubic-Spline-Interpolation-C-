/**
 * @file    Matrix.h
 * @date    Created on 10/7/2020 at 17:01
 * @author  Ali Karimiafshar <Karimiafsharali@gmail.com>
 * @brief   
 */

#ifndef MATRIX_MULTIPLIER_MATRIX_H
#define MATRIX_MULTIPLIER_MATRIX_H

#include "General_Include.h"

template <class Object>
class Matrix {
  protected:
	int rows;
	int cols;
	Object** grid;
  
  public:
	explicit Matrix(int row = 1, int col = 1);
	
	virtual ~Matrix();
	
	int getRow() const;
	
	int getCol() const;
	
	Object** getGrid() const;

	Object getElement(int row, int col);

	void setElement(int row, int col, Object const &value);
	
	void printMatrix();
	
	void setValFromUser();
	
	void setVal(int row, int col, const Object &val);
	
	Matrix mMultiply(const Matrix<Object> &M);
	
	void setToZero();
	
	void setTridiagonal();
	
	Object* operator[] (int index);
	
}; // End of class Matrix.h

#endif //MATRIX_MULTIPLIER_MATRIX_H