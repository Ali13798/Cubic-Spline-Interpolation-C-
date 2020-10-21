/**
 * @file    Matrix.cpp
 * @date    Created on 10/7/2020 at 17:01
 * @author  Ali Karimiafshar <Karimiafsharali@gmail.com>
 * @brief   
 */

#include "Matrix.h"

template
class Matrix<int>;

template
class Matrix<float>;

template <class Object>
Matrix<Object>::Matrix(int row, int col) {     // Constructor
	rows = row;
	cols = col;
	
	// Create the rows
	grid = new Object* [rows];
	// Create the columns
	for (int i{0}; i < rows; ++i) {
		grid[i] = new Object[cols];
	}
}

template <class Object>
Matrix<Object>::~Matrix() {
	for (int i = 0; i < rows; ++i) {
		delete[] grid[i];
	}
	delete[] grid;
}

template <class Object>
int Matrix<Object>::getRow() const {
	return rows;
}

template <class Object>
int Matrix<Object>::getCol() const {
	return cols;
}

template <class Object>
Object** Matrix<Object>::getGrid() const {
	return grid;
}

template <class Object>
Object Matrix<Object>::getElement(int row, int col) {
	return grid[row][col];
}

template <class Object>
void Matrix<Object>::setElement(int row, int col, Object const &value) {
	grid[row][col] = value;
}

template <class Object>
void Matrix<Object>::printMatrix() {
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			cout << grid[i][j] << " ";
		}
		cout << '\n';
	}
	cout << '\n';
}

template <class Object>
void Matrix<Object>::setValFromUser() {
	for (int row = 0; row < rows; row++) {
		for (int col = 0; col < cols; col++) {
			Object tempVal;
			cin >> tempVal;
			grid[row][col] = tempVal;
		}
	}
}

template <class Object>
void Matrix<Object>::setVal(int row, int col, const Object &val) {
	grid[row][col] = val;
}

template <class Object>
Matrix<Object> Matrix<Object>::mMultiply(const Matrix<Object> &M) {
	Matrix<Object> r(rows, M.getCol());
	
	if (cols != M.getRow()) {
		cout << "\nCannot multiply these matrices...";
		exit(EXIT_FAILURE);
	} else {
		r.setToZero();
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < M.getCol(); j++) {
				for (int k = 0; k < cols; k++) {
					Object tempVal = r.getGrid()[i][j];
					tempVal += Matrix::grid[i][k]*M.getGrid()[k][j];
					r.setVal(i, j, tempVal);
				}
			}
		}
		return r;
	}
}

template <class Object>
void Matrix<Object>::setToZero() {
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			Matrix::grid[i][j] = 0;
		}
	}
}

template <class Object>
void Matrix<Object>::setTridiagonal() {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (i == j) {
				grid[i][j] = 4;
			} else if (i > 0 && j == i - 1) {
				grid[i][j] = 1;
			} else if (j > 0 && i == j - 1) {
				grid[i][j] = 1;
			} else {
				grid[i][j] = 0;
			}
		}
	}
}

template <class Object>
Object* Matrix<Object>::operator[](int index) {
	return grid[index];
}