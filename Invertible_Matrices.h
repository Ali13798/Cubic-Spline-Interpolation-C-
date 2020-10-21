/**
 * @file    Invertible_Matrices.h
 * @date    Created on 10/12/2020 at 12:58
 * @author  Ali Karimiafshar <Karimiafsharali@gmail.com>
 * @brief   
 */

#ifndef MATRIX_MULTIPLIER_INVERTIBLE_MATRICES_H
#define MATRIX_MULTIPLIER_INVERTIBLE_MATRICES_H

#include "Matrix.h"


class InvertibleMatrix : public Matrix<int> {
  private:
	const int MSIZE;
  public:
	explicit InvertibleMatrix(const int msize);
	virtual ~InvertibleMatrix() = default;

	/**
 	 * The code was copied form the internet to find the inverse of a matrix
 	 * https://www.geeksforgeeks.org/adjoint-inverse-matrix/
 	 */

	// Function to get cofactor of A[p][q] in temp[][]. n is current dimension of A[][]
	void getCofactor(InvertibleMatrix &temp, int p, int q, int n);

	// Recursive function for finding determinant of matrix. n is current dimension of A[][].
	int determinant(int n);

	// Function to get adjoint of A[N][N] in adj[N][N].
	void adjoint(Matrix<int> &adj);

	// Function to calculate and store inverse, returns false if matrix is singular
	bool inverse(Matrix<float> &inverse);


}; // End of class Invertible_Matrices.h




#endif //MATRIX_MULTIPLIER_INVERTIBLE_MATRICES_H
