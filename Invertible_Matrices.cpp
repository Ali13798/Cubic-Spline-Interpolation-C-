/**
 * @file    Invertible_Matrices.cpp
 * @date    Created on 10/12/2020 at 12:58
 * @author  Ali Karimiafshar <Karimiafsharali@gmail.com>
 * @brief   
 */

#include "Invertible_Matrices.h"

InvertibleMatrix::InvertibleMatrix(const int msize) : MSIZE(rows), Matrix(msize, msize) {
	// Intentionally Empty
}


/**
 * The code was copied form the internet to find the inverse of a matrix
 * https://www.geeksforgeeks.org/adjoint-inverse-matrix/
 */

// Function to get cofactor of A[p][q] in temp[][]. n is current dimension of A[][]
void InvertibleMatrix::getCofactor(InvertibleMatrix &temp, int p, int q, int n) {
	int i = 0, j = 0;

	// Looping for each element of the matrix
	for (int row = 0; row < n; row++) {
		for (int col = 0; col < n; col++) {
			//  Copying into temporary matrix only those element
			//  which are not in given row and column
			if (row != p && col != q) {
				temp[i][j++] = grid[row][col];

				// Row is filled, so increase row index and
				// reset col index
				if (j == n - 1) {
					j = 0;
					i++;
				}
			}
		}
	}
}

// Recursive function for finding determinant of matrix. n is current dimension of A[][].
int InvertibleMatrix::determinant(int n) {
	int D = 0; // Initialize result

	//  Base case : if matrix contains single element
	if (n == 1) {
		return grid[0][0];
	}

	InvertibleMatrix *temp = new InvertibleMatrix(MSIZE); 
	// InvertibleMatrix temp(MSIZE); // To store cofactors

	int sign = 1;  // To store sign multiplier

	// Iterate for each element of first row
	for (int f = 0; f < n; f++) {
		// Getting Cofactor of A[0][f]
		getCofactor(*temp, 0, f, n);
		D += sign*grid[0][f]*temp->determinant(n - 1);

		// terms are to be added with alternate sign
		sign = -sign;
	}

	return D;
}

// Function to get adjoint of A[MSIZE][MSIZE] in adj[MSIZE][MSIZE].
void InvertibleMatrix::adjoint(Matrix<int> &adj) {
	if (MSIZE == 1) {
		adj[0][0] = 1;
		return;
	}

	// temp is used to store cofactors of A[][]
	int sign = 1;
	InvertibleMatrix *temp = new InvertibleMatrix(MSIZE);
	// InvertibleMatrix temp(MSIZE);

	for (int i = 0; i < MSIZE; i++) {
		for (int j = 0; j < MSIZE; j++) {
			// Get cofactor of A[i][j]
			getCofactor(*temp, i, j, MSIZE);

			// sign of adj[j][i] positive if sum of row
			// and column indexes is even.
			sign = ( ( i + j )%2 == 0 ) ? 1 : -1;

			// Interchanging rows and columns to get the
			// transpose of the cofactor matrix
			adj[j][i] = ( sign )*( temp->determinant(MSIZE - 1) );
		}
	}
}

// Function to calculate and store inverse, returns false if matrix is singular
bool InvertibleMatrix::inverse(Matrix<float> &inverse) {
	// Find determinant of A[][]
	int det = determinant(MSIZE);
	if (det == 0) {
		cout << "Singular matrix, can't find its inverse";
		return false;
	}

	// Find adjoint
	InvertibleMatrix *adj = new InvertibleMatrix(MSIZE);
	// InvertibleMatrix adj(MSIZE);
	adjoint(*adj);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
	for (int i = 0; i < MSIZE; i++) {
		for (int j = 0; j < MSIZE; j++) {
			inverse[i][j] = adj->getElement(i, j)/float(det);
		}
	}

	return true;
}