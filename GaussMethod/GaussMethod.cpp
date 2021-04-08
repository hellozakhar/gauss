#include <iostream>

#include "Utils.h"
#include "Timer.h"
#include "GaussMethod.h"

const double EPS = 1e-5;

void gaussDef(double* a, double* y, double* x, int n) {
	double mu = 0.0;

	for (int i = 0; i < n - 1; i++) {
		for (int k = i + 1; k <= n - 1; k++) {

			mu = a[k * n + i] / a[i * n + i];

			for (int j = i; j <= n - 1; j++) {
				a[k * n + j] = a[k * n + j] - mu * a[i * n + j];

			}

			y[k] = y[k] - mu * y[i];
		}
	}

	for (int i = n - 1; i >= 0; i--) {
		x[i] = y[i];

		for (int j = i + 1; j < n; j++) {
			x[i] -= a[i * n + j] * x[j];
		}
		
		x[i] /= a[i * n + i];
	}
}

void gaussMax(double* a, double* y, double* x, int n) {
	double max = 0.0;
	double mu = 0.0;
	int index = 0;
	int* row = new int[n];

	for (int i = 0; i < n; i++) {
		row[i] = i;
	}

	for (int i = 0; i < n - 1; i++) {

		max = a[i * n + i];
		index = i;

		for (int j = i + 1; j < n; j++) {
			if (abs(a[j * n + i]) > max) {
				max = abs(a[j * n + i]);
				index = j;
			}
		}

		swap(row[i], row[index]);

		for (int k = i + 1; k <= n - 1; k++) {

			mu = a[row[k] * n + i] / a[row[i] * n + i];

			for (int j = i; j <= n - 1; j++) {
				a[row[k] * n + j] = a[row[k] * n + j] - mu * a[row[i] * n + j];

			}

			y[row[k]] = y[row[k]] - mu * y[row[i]];
		}
	}

	for (int i = n - 1; i >= 0; i--) {
		x[i] = y[row[i]];

		for (int j = i + 1; j < n; j++) {
			x[i] -= a[row[i] * n + j] * x[j];
		}

		x[i] /= a[row[i] * n + i];
	}
}

void LU_Decomposition(double* A, double* L, double* U, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			U[i * n + j] = A[i * n + j];
		}
	}

	for (int i = 0; i < n; i++) {
		L[i * n + i] = 1;

		for (int k = i + 1; k < n; k++) {
			double mu = U[k * n + i] / U[i * n + i];

			for (int j = i; j < n; j++) {
				U[k * n + j] -= mu * U[i * n + j];
			}

			L[k * n + i] = mu;
			L[i * n + k] = 0;
		}
	}

	for (int i = 1; i < n; i++) {
		for (int j = 0; j < i; j++) {
			U[i * n + j] = 0;
		}
	}
}

void LU_Solve(double* L, double* U, double* y, double* x, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			y[i] -= L[i * n + j] * y[j];
		}
	}

	for (int i = n - 1; i >= 0; i--) {
		x[i] = y[i];

		for (int j = i + 1; j < n; j++) {
			x[i] -= U[i * n + j] * x[j];
		}

		x[i] /= U[i * n + i];
	}
}

void LU_Decomposition(double* A, int n) {
	for (int i = 0; i < n; i++) {
		for (int k = i + 1; k < n; k++) {
			double mu = A[k * n + i] / A[i * n + i];

			for (int j = i; j < n; j++) {
				A[k * n + j] -= mu * A[i * n + j];
			}

			A[k * n + i] = mu;
		}
	}
}

void LU_Solve(double* LU, double* y, double* x, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			y[i] -= LU[i * n + j] * y[j];
		}
	}

	for (int i = n - 1; i >= 0; i--) {
		x[i] = y[i];

		for (int j = i + 1; j < n; j++) {
			x[i] -= LU[i * n + j] * x[j];
		}

		x[i] /= LU[i * n + i];
	}
}

void testGaussDef(bool isAuto) {
	int n = 0;
	std::cout << "Enter count of equations: ";
	std::cin >> n;

	double* a = new double[n * n];
	double* y = new double[n];
	double* x = new double[n];
	init(x, n);

	if (isAuto) {
		generateSystemOfLinearEquations(a, y, n);
	}
	else {
		input(a, y, n);
	}

	//print(a, y, n);
	gaussDef(a, y, x, n);

	for (int i = 0; i < n; i++) {
		std::cout << "x[" << i << "]=" << x[i] << std::endl;
	}

	delete[] a;
	delete[] y;
	delete[] x;
}

void testGaussMax(bool isAuto) {
	int n;
	std::cout << "Enter count of equations: ";
	std::cin >> n;

	double* a = new double[n * n];
	double* y = new double[n];
	double* x = new double[n];
	init(x, n);

	if (isAuto) {
		generateSystemOfLinearEquations(a, y, n);
	}
	else {
		input(a, y, n);
	}

	gaussMax(a, y, x, n);

	for (int i = 0; i < n; i++) {
		std::cout << "x[" << i << "]=" << x[i] << std::endl;
	}

	delete[] a;
	delete[] y;
	delete[] x;
}

void testGaussDiscrepancy() {
	int n = 0;
	std::cout << "Enter count of equations: ";
	std::cin >> n;

	double* aDef = new double[n * n];
	double* yDef = new double[n];
	double* xDef = new double[n];
	init(xDef, n);
	double* aMax = new double[n * n];
	double* yMax = new double[n];
	double* xMax = new double[n];
	init(xMax, n);
	double* aOrig = new double[n * n];
	double* yOrig = new double[n];

	generateSystemOfLinearEquations(aDef, yDef, n);

	copy(aDef, aMax, n * n);
	copy(yDef, yMax, n * n);
	copy(aDef, aOrig, n * n);
	copy(yDef, yOrig, n * n);

	gaussDef(aDef, yDef, xDef, n);
	gaussMax(aMax, yMax, xMax, n);

	delete[] aDef;
	delete[] yDef;
	delete[] aMax;
	delete[] yMax;

	for (int i = 0; i < n; i++) {
		std::cout << "(def) x[" << i << "]=" << xDef[i] << std::endl;
	}

	for (int i = 0; i < n; i++) {
		std::cout << "(max) x[" << i << "]=" << xMax[i] << std::endl;
	}

	std::cout << "\nDiscrepancy" << std::endl;
	double* discrepancyGaussDef = new double[n];
	double* discrepancyGaussMax = new double[n];

	for (int i = 0; i < n; i++) {
		discrepancyGaussDef[i] = yOrig[i];
		discrepancyGaussMax[i] = yOrig[i];

		for (int j = 0; j < n; j++) {
			discrepancyGaussDef[i] -= aOrig[i * n + j] * xDef[j];
			discrepancyGaussMax[i] -= aOrig[i * n + j] * xMax[j];
		}
	}

	for (int i = 0; i < n; i++) {
		std::cout << "(discrepancy def) delta[" << i << "]=" << discrepancyGaussDef[i] << std::endl;
		std::cout << "(discrepancy max) delta[" << i << "]=" << discrepancyGaussMax[i] << std::endl;
	}
	delete[] xDef;
	delete[] xMax;

	delete[] aOrig;
	delete[] yOrig;

	delete[] discrepancyGaussDef;
	delete[] discrepancyGaussMax;
}

void testLUDecomposition() {
	int n = 0;
	std::cout << "Enter count of equations: ";
	std::cin >> n;

	double* A = new double[n * n];
	double* L = new double[n * n];
	double* U = new double[n * n];
	init(A, n * n);
	init(L, n * n);
	init(U, n * n);

	generateRandomMatrix(A, n);
	printMatrix(A, n);
	LU_Decomposition(A, L, U, n);
	printMatrix(L, n);
	printMatrix(U, n);

	delete[] A;
	delete[] L;
	delete[] U;
}

void testCompareSolutionGaussAndLU() {
	int n = 3;
	double* A_GAUSS = new double[n * n];
	double* Y_GAUSS = new double[n];
	double* X_GAUSS = new double[n];

	double* A_GAUSS_MAX = new double[n * n];
	double* Y_GAUSS_MAX = new double[n];
	double* X_GAUSS_MAX = new double[n];

	double* A_OLDLU = new double[n * n];
	double* L_OLDLU = new double[n * n];
	double* U_OLDLU = new double[n * n];
	double* Y_OLDLU = new double[n];
	double* X_OLDLU = new double[n];

	double* A_LU = new double[n * n];
	double* Y_LU = new double[n];
	double* X_LU = new double[n];

	/*generateSystemOfLinearEquations(A_GAUSS, Y_GAUSS, n);*/
	A_GAUSS[0] = 2; A_GAUSS[1] = 1; A_GAUSS[2] = 4;
	A_GAUSS[3] = 3; A_GAUSS[4] = 2; A_GAUSS[5] = 1;
	A_GAUSS[6] = 1; A_GAUSS[7] = 3; A_GAUSS[8] = 3;
	Y_GAUSS[0] = 16; Y_GAUSS[1] = 10; Y_GAUSS[2] = 16;

	copy(A_GAUSS, A_LU, n * n);
	copy(A_GAUSS, A_OLDLU, n * n);
	copy(A_GAUSS, A_GAUSS_MAX, n * n);
	copy(Y_GAUSS, Y_LU, n);
	copy(Y_GAUSS, Y_OLDLU, n);
	copy(Y_GAUSS, Y_GAUSS_MAX, n);

	printMatrix(A_GAUSS, n);
	print(Y_GAUSS, n);

	gaussDef(A_GAUSS, Y_GAUSS, X_GAUSS, n);
	gaussMax(A_GAUSS_MAX, Y_GAUSS_MAX, X_GAUSS_MAX, n);
	LU_Decomposition(A_OLDLU, L_OLDLU, U_OLDLU, n);
	LU_Decomposition(A_LU, n);
	LU_Solve(L_OLDLU, U_OLDLU, Y_OLDLU, X_OLDLU, n);
	LU_Solve(A_LU, Y_LU, X_LU, n);

	print(X_GAUSS, n);
	print(X_GAUSS_MAX, n);
	print(X_OLDLU, n);
	print(X_LU, n);

	delete[] A_GAUSS;
	delete[] Y_GAUSS;
	delete[] X_GAUSS;
	delete[] A_GAUSS_MAX;
	delete[] Y_GAUSS_MAX;
	delete[] X_GAUSS_MAX;
	delete[] A_LU;
	delete[] Y_LU;
	delete[] X_LU;
	delete[] A_OLDLU;
	delete[] L_OLDLU;
	delete[] U_OLDLU;
	delete[] Y_OLDLU;
	delete[] X_OLDLU;

}

void test_OLDLU_solve(int n) {
	double* A = new double[n * n];
	double* Y = new double[n];
	double* L = new double[n * n];
	double* U = new double[n * n];
	double* X = new double[n];

	generateSystemOfLinearEquations(A, Y, n);
	LU_Decomposition(A, L, U, n);
	LU_Solve(L, U, Y, X, n);

	delete[] A;
	delete[] Y;
	delete[] L;
	delete[] U;
	delete[] X;
}

void test_LU_solve(int n) {
	double* A = new double[n * n];
	double* Y = new double[n];
	double* X = new double[n];

	generateSystemOfLinearEquations(A, Y, n);
	LU_Decomposition(A, n);
	LU_Solve(A, Y, X, n);

	delete[] A;
	delete[] Y;
	delete[] X;
}

int main() {
	// isAuto true  => coeffs and y - generated,
	// isAuto false => coeffs and y - enter manually.
	const bool isAuto = true;
	const int n = 2500;

	//testGaussDef(isAuto);
	//testGaussDiscrepancy();
	//testLUDecomposition();
	//testCompareSolutionGaussAndLU();

	//test_OLDLU_solve(n);
	//test_LU_solve(n);


	int zzz;
	std::cin >> zzz;
	return 0;
}