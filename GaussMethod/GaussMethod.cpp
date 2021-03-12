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

	copy(aDef, aMax, n, n);
	copy(yDef, yMax, n, n);
	copy(aDef, aOrig, n, n);
	copy(yDef, yOrig, n, n);

	gaussDef(aDef, yDef, xDef, n);
	gaussMax(aMax, yMax, xMax, n);

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

	delete[] aDef;
	delete[] yDef;
	delete[] xDef;

	delete[] aMax;
	delete[] yMax;
	delete[] xMax;

	delete[] aOrig;
	delete[] yOrig;

	delete[] discrepancyGaussDef;
	delete[] discrepancyGaussMax;
}

int main() {
	// isAuto true  => coeffs and y - generated,
	// isAuto false => coeffs and y - enter manually.
	const bool isAuto = true;

	//testGaussDef(isAuto);
	testGaussDiscrepancy();

	return 0;
}