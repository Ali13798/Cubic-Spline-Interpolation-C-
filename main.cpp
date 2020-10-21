/**
 * @file    main.cpp
 * @date    Created on 10/12/2020 at 12:17
 * @author  Ali Karimiafshar <Karimiafsharali@gmail.com>
 * @brief	Calculate the cubic spline interpolation
 */

#include "General_Include.h"
#include "Matrix.h"
#include "Invertible_Matrices.h"

float calc_Y(float dX, Matrix<float>* coefficients, int posCounter);

Matrix<float>* calcInterpolate(int numPoints, int subIntervalPoints, float xValues[], Matrix<float>* coefficients);

int main() {
	int numPoints;
	cout << "Enter the number of points desired: ";
	cin >> numPoints;
	
	cout << "Enter the point coordinates in order, with each x and y value separated ";
	cout << "with a space and each point coordinate separated with a space:\n";
	float xValues[numPoints];
	float yValues[numPoints];
	for (int i = 0; i < numPoints; i++) {
		cin >> xValues[i] >> yValues[i];
	}
	
	cout << "\nThe points entered are:\n";
	for (int i = 0; i < numPoints; i++) {
		cout << xValues[i] << ", " << yValues[i] << '\n';
	}
	
	cout << "\nThe tridiagonal matrix is:\n";
	int mSize = numPoints - 2;
	auto* tridiagonal = new InvertibleMatrix(mSize);
	tridiagonal->setTridiagonal();
	tridiagonal->printMatrix();
	
	cout << "The inverse tridiagonal matrix is: \n";
	auto* inverseTridiagonal = new Matrix<float>(mSize, mSize);
	tridiagonal->inverse(*inverseTridiagonal);
	inverseTridiagonal->printMatrix();
	delete tridiagonal;
	
	cout << "The Y equation matrix is: \n";
	float hValue = xValues[1] - xValues[0];
	auto* yMatrix = new Matrix<float>(mSize, 1);
	for (int i = 0; i < mSize; i++) {
		float val = 6/( hValue*hValue )*( yValues[i] - 2*yValues[i + 1] + yValues[i + 2] );
		yMatrix->setElement(i, 0, val);
	}
	yMatrix->printMatrix();
	
	cout << "The M matrix is: \n";
	auto theMMatrix = inverseTridiagonal->mMultiply(*yMatrix);
	theMMatrix.printMatrix();
	
	//// A, B, C, D, mValues
	auto* coefficients = new Matrix<float>(5, numPoints);
	
	coefficients->setElement(4, 0, 0);
	coefficients->setElement(4, numPoints - 1, 0);
	for (int i = 1; i < numPoints - 1; i++) {
		coefficients->setElement(4, i, theMMatrix.getElement(i - 1, 0));
	}
	
	for (int i = 0; i < numPoints; i++) {
		float A_coeff = ( coefficients->getElement(4, i + 1) - coefficients->getElement(4, i) )/( 6*hValue );
		coefficients->setElement(0, i, A_coeff);
		
		float B_Coeff = ( coefficients->getElement(4, i)/2 );
		coefficients->setElement(1, i, B_Coeff);
		
		float C_coeff = ( yValues[i + 1] - yValues[i] )/hValue;
		C_coeff -= hValue*( coefficients->getElement(4, i + 1) + 2*coefficients->getElement(4, i) )/6;
		coefficients->setElement(2, i, C_coeff);
		
		float D_Coeff = yValues[i];
		coefficients->setElement(3, i, D_Coeff);
	}
	
	int menuChoice;
	cout << "Select from the menu:\n";
	cout << "\t1. Enter the number of points in each interval.\n";
	cout << "\t2. Enter the maximum error desired.\n";
	cout << "Enter choice: ";
	cin >> menuChoice;
	
	int subIntervalPoints = 10;
	if (menuChoice == 1) {
		cout << "Enter How many points in each interval (excluding boundaries): ";
		cin >> subIntervalPoints;
		int numFinalCoordinates = subIntervalPoints*( numPoints - 1 ) + numPoints;
		
		auto interpolatedPoints = calcInterpolate(numPoints, subIntervalPoints, xValues, coefficients);
		
		cout << "The Cubic Spline Interpolation calculation is done.\n";
		
		ofstream myFile;
		myFile.open("Ali_CSI_Points.txt", std::ios::out);
		for (int i = 0; i < numFinalCoordinates; i++) {
			myFile << interpolatedPoints->getElement(0, i) << " " << interpolatedPoints->getElement(1, i) << '\n';
		}
		myFile.close();
		
		ofstream AML_File;
		AML_File.open("Ali_CSI.aml", std::ios::out);
		for (int i = 0; i < numFinalCoordinates; i++) {
			AML_File << "PMOVE(" << interpolatedPoints->getElement(0, i) << ", "
			         << interpolatedPoints->getElement(1, i)
			         << ", 0, 0);\n";
		}
		AML_File.close();
		
	} else if (menuChoice == 2) {
		float userError;
		cout << "Enter the max error acceptable: ";
		cin >> userError;
		
		int numFinalCoordinates;
		int loopEnd = 0;
		while (loopEnd < 10) {
			numFinalCoordinates = subIntervalPoints*( numPoints - 1 ) + numPoints;
			auto* interpolated_X = new float[numFinalCoordinates];
			auto* interpolated_Y = new float[numFinalCoordinates];
			float delta = hValue/( subIntervalPoints + 1 );
			
			float inputX = xValues[0];
			int posCounter = 0;
			int subInterval = 0;
			
			for (int i = 0; i < numFinalCoordinates; i++) {
				float X_minus_Xi = inputX - xValues[posCounter];
				interpolated_Y[i] = calc_Y(X_minus_Xi, coefficients, posCounter);
				interpolated_X[i] = inputX;
				
				inputX += delta;
				subInterval++;
				
				if (subInterval > subIntervalPoints) {
					posCounter++;
					inputX = xValues[posCounter];
					subInterval = 0;
				}
			}
			
			loopEnd++;
		}
		
		cout << "The Cubic Spline Interpolation calculation is done.\n";

//		ofstream myFile;
//		myFile.open("Ali_CSI_Points.txt", std::ios::out);
//		for (int i = 0; i < numFinalCoordinates; i++) {
//			myFile << interpolated_X[i] << " " << interpolated_Y[i] << '\n';
//		}
//		myFile.close();
//
//		ofstream AML_File;
//		AML_File.open("Ali_CSI.aml", std::ios::out);
//		for (int i = 0; i < numFinalCoordinates; i++) {
//			AML_File << "PMOVE(" << interpolated_X[i] << ", " << interpolated_Y[i];
//			AML_File << ", 0, 0);\n";
//		}
//		AML_File.close();
	}
	
	return EXIT_SUCCESS;
}

float calc_Y(float dX, Matrix<float>* coefficients, int posCounter) {
	float y = 0;
	y += coefficients->getElement(0, posCounter)*dX*dX*dX;
	y += coefficients->getElement(1, posCounter)*dX*dX;
	y += coefficients->getElement(2, posCounter)*dX;
	y += coefficients->getElement(3, posCounter);
	return y;
}

Matrix<float>* calcInterpolate(int numPoints, int subIntervalPoints, float xValues[], Matrix<float>* coefficients) {
	int numFinalCoordinates = subIntervalPoints*( numPoints - 1 ) + numPoints;
	auto M = new Matrix<float>(2, numFinalCoordinates);
	
	float hValue = xValues[1] - xValues[0];
	float delta = hValue/( subIntervalPoints + 1 );
	
	float inputX = xValues[0];
	int posCounter = 0;
	int subInterval = 0;
	
	for (int i = 0; i < numFinalCoordinates; i++) {
		float X_minus_Xi = inputX - xValues[posCounter];
		auto y = calc_Y(X_minus_Xi, coefficients, posCounter);
		
		M->setElement(0, i, inputX);
		M->setElement(1, i, y);
		
		inputX += delta;
		subInterval++;
		
		if (subInterval > subIntervalPoints) {
			posCounter++;
			inputX = xValues[posCounter];
			subInterval = 0;
		}
	}
	return M;
}