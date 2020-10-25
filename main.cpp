/**
 * @file    main.cpp
 * @date    Created on 10/12/2020 at 12:17
 * @author  Ali Karimiafshar <Karimiafsharali@gmail.com>
 * @brief	Calculate the cubic spline interpolation
 */

#include "General_Include.h"
#include "Matrix.h"
#include "Invertible_Matrices.h"

// prototypes:

Matrix<float>* calcInterpolate(int numPts, int subPts, const float xValues[], Matrix<float>* coefficients);

float calc_Y(float dX, Matrix<float>* coefficients, int posCounter);

void writeOutput(int points, Matrix<float>* M);


int main() {
	int numPts;
	cout << "Enter the number of points desired: ";
	cin >> numPts;
	
	cout << "Enter the point coordinates in order, with each x and y value separated "
	     << "with a space and each point coordinate separated with a space:\n";
	float xValues[numPts];
	float yValues[numPts];
	for (int i = 0; i < numPts; i++) {
		cin >> xValues[i] >> yValues[i];
	}
	
	cout << "\nThe points entered are:\n";
	for (int i = 0; i < numPts; i++) {
		cout << xValues[i] << ", " << yValues[i] << '\n';
	}
	
	cout << "\nThe tridiagonal matrix is:\n";
	int mSize = numPts - 2;
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
	auto* coefficients = new Matrix<float>(5, numPts);
	
	coefficients->setElement(4, 0, 0);
	coefficients->setElement(4, numPts - 1, 0);
	for (int i = 1; i < numPts - 1; i++) {
		coefficients->setElement(4, i, theMMatrix.getElement(i - 1, 0));
	}
	
	for (int i = 0; i < numPts; i++) {
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
	cout << "\t3. Exit the program.\n";
	cout << "Enter choice: ";
	cin >> menuChoice;
	
	int subPts = 1;
	if (menuChoice == 1) {
		cout << "Enter How many points in each interval (excluding boundaries): ";
		cin >> subPts;
		int numFinalPts = ( numPts - 1 )*( subPts + 1 ) + 1;
		
		auto interpolatedPts = calcInterpolate(numPts, subPts, xValues, coefficients);
		
		writeOutput(numFinalPts, interpolatedPts);
		
		cout << "The Cubic Spline Interpolation calculation is done.\n";
		
		return EXIT_SUCCESS;
		
	} else if (menuChoice == 2) {
		float userError;
		cout << "Enter the max error acceptable: ";
		cin >> userError;
		
		bool done = false;
		while (!done) {
			int subSubPts = 100; // how many intervals is each subinterval divided into to calculate the error
			int subPtsEq = ( subSubPts + 1 )*( subPts + 1 ) - 1; // Subinterval points equivalent
			auto tempPoints = calcInterpolate(numPts, subPts, xValues, coefficients);
			auto truePoints = calcInterpolate(numPts, subPtsEq, xValues, coefficients);
			
			int numFinalPts = ( numPts - 1 )*( subPtsEq + 1 ) + 1;
			auto result = new Matrix<float>(3, numFinalPts);
			
			float true_error = 0;
			for (int interval = 0, xPos = 0, slopePos = 0, posCounter = 0; interval < numPts - 1; interval++) {
				// for each interval
				
				for (int subInt = 0; subInt < subPts + 1; subInt++, slopePos++, xPos++) {
					// for each subinterval
					auto x1 = tempPoints->getElement(0, slopePos);
					auto x2 = tempPoints->getElement(0, slopePos + 1);
					auto y1 = tempPoints->getElement(1, slopePos);
					auto y2 = tempPoints->getElement(1, slopePos + 1);
					auto slope = ( y2 - y1 )/( x2 - x1 );
					
					float x_val = tempPoints->getElement(0, xPos);
					float y_val = tempPoints->getElement(1, xPos);
					
					for (int subSubInt = 0; subSubInt < subSubPts + 1; subSubInt++, posCounter++) {
						// for each sub subinterval
						float delta = hValue/( subPtsEq + 1 );
						
						float true_x = truePoints->getElement(0, posCounter);
						float true_y = truePoints->getElement(1, posCounter);
						float line_y = y_val + subSubInt*delta*slope;
						
						result->setElement(0, posCounter, true_x);
						result->setElement(1, posCounter, line_y);
						result->setElement(2, posCounter, true_y);
						
						float error = std::abs(line_y - true_y);
						true_error = std::max(error, true_error);
						
					}
				}
			}
			
			if (true_error < userError) {
				done = true;
				cout << "SUBPOINTS REQ: " << subPts << endl
				     << "MAX ERROR: " << true_error << endl;
				
				result->setElement(0, numFinalPts - 1, xValues[numPts - 1]);
				result->setElement(1, numFinalPts - 1, yValues[numPts - 1]);
				result->setElement(2, numFinalPts - 1, yValues[numPts - 1]);
				
				writeOutput(numFinalPts, result);
			}
			subPts += 2;
		}
		
		cout << "The Cubic Spline Interpolation calculation is done.\n";
		
		return EXIT_SUCCESS;
	} else {
		cout << "Goodbye.\n";
		return EXIT_SUCCESS;
	}
}

Matrix<float>*
calcInterpolate(int numPts, int subPts, const float xValues[], Matrix<float>* coefficients) {
	int numFinalPts = ( numPts - 1 )*( subPts + 1 ) + 1;
	auto M = new Matrix<float>(2, numFinalPts);
	
	float hValue = xValues[1] - xValues[0];
	float delta = hValue/( subPts + 1 );
	
	float inputX = xValues[0];
	int posCounter = 0;
	int subInterval = 0;
	
	for (int i = 0; i < numFinalPts; i++) {
		float X_minus_Xi = inputX - xValues[posCounter];
		auto y = calc_Y(X_minus_Xi, coefficients, posCounter);
		
		M->setElement(0, i, inputX);
		M->setElement(1, i, y);
		
		inputX += delta;
		subInterval++;
		
		if (subInterval > subPts) {
			posCounter++;
			inputX = xValues[posCounter];
			subInterval = 0;
		}
	}
	return M;
}

float calc_Y(float dX, Matrix<float>* coefficients, int posCounter) {
	float y = 0;
	y += coefficients->getElement(0, posCounter)*dX*dX*dX;
	y += coefficients->getElement(1, posCounter)*dX*dX;
	y += coefficients->getElement(2, posCounter)*dX;
	y += coefficients->getElement(3, posCounter);
	return y;
}

void writeOutput(int points, Matrix<float>* M) {
	ofstream myFile;
	myFile.open("Ali_CSI_Points.txt", std::ios::out);
	for (int i = 0; i < points; i++) {
		myFile << M->getElement(0, i) << " " << M->getElement(1, i) << '\n';
	}
	myFile.close();
	
	ofstream AML_File;
	AML_File.open("Ali_CSI.aml", std::ios::out);
	for (int i = 0; i < points; i++) {
		AML_File << "PMOVE(PT(" << M->getElement(0, i) << ", "
		         << M->getElement(1, i)
		         << ", 0, 0));" << endl;
	}
	AML_File.close();
}