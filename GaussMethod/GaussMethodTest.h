#pragma once
#include <iostream>

#include "Timer.h"
#include "GaussMethodAux.h"
#include "GaussMethod.h"

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
	decompose_v1(A_OLDLU, L_OLDLU, U_OLDLU, n);
	decompose_v2_par(A_LU, n);
	solve(L_OLDLU, U_OLDLU, Y_OLDLU, X_OLDLU, n);
	solve(A_LU, Y_LU, X_LU, n);

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

void testDecomposeSolveV1(int n) {
	double* A = new double[n * n];
	double* Y = new double[n];
	double* L = new double[n * n];
	double* U = new double[n * n];
	double* X = new double[n];

	generateSystemOfLinearEquations(A, Y, n);
	decompose_v1(A, L, U, n);
	solve(L, U, Y, X, n);
	std::cin >> n;
	delete[] A;
	delete[] Y;
	delete[] L;
	delete[] U;
	delete[] X;
}

void testDecomposeSolveV2(int n) {
	double* A = new double[n * n];
	double* Y = new double[n];
	double* X = new double[n];

	generateSystemOfLinearEquations(A, Y, n);
	decompose_v2(A, n);
	solve(A, Y, X, n);
	std::cin >> n;
	delete[] A;
	delete[] Y;
	delete[] X;
}

void testDecomposeSolveV1(int n, int times, bool isForCopy) {
	if (!isForCopy) {
		std::cout << "\nLU-decomposition (2ARR) of matrix  " << n << " * " << n
			<< "\nRun " << times << " times." << std::endl;
	}

	double* A = new double[n * n];
	double* Y = new double[n];
	double* L = new double[n * n];
	double* U = new double[n * n];
	double* X = new double[n];
	double* time = new double[times];

	generateSystemOfLinearEquations(A, Y, n);

	for (int i = 0; i < times; i++) {

		/*init(L, n*n);
		init(U, n*n);
		init(Y, n*n);
		init(X, n*n);*/

		Timer t1;
		decompose_v1(A, L, U, n);
		solve(L, U, Y, X, n);
		time[i] = t1.elapsed();
	}

	double avg = 0.0;
	double min = time[0];
	double max = time[0];

	for (int i = 0; i < times; i++) {
		//std::cout << i << "\ttime taken: " << time.at(i) << '\n';
		avg += time[i];
		if (min > time[i]) min = time[i];
		if (max < time[i]) max = time[i];
	}
	avg /= times;

	if (!isForCopy) {
		std::cout << "------TIME------"
			<< "\nAVG: " << avg
			<< "\nMIN: " << min
			<< "\nMAX: " << max << std::endl;
	}
	else {
		std::cout << avg << std::endl;
	}

	delete[] A;
	delete[] Y;
	delete[] L;
	delete[] U;
	delete[] X;
	delete[] time;
}

void testDecomposeSolveV2(int n, int times, bool isForCopy) {
	if (!isForCopy) {
		std::cout << "\nLU-decomposition (1ARR) of matrix  " << n << " * " << n
			<< "\nRun " << times << " times." << std::endl;
	}

	double* A = new double[n * n];
	double* Y = new double[n];
	double* X = new double[n];
	double* time = new double[times];


	for (int i = 0; i < times; i++) {
		generateSystemOfLinearEquations(A, Y, n);
		/*init(L, n*n);
		init(U, n*n);
		init(Y, n*n);
		init(X, n*n);*/

		Timer t1;
		decompose_v2(A, n);
		solve(A, Y, X, n);
		time[i] = t1.elapsed();
	}

	double avg = 0.0;
	double min = time[0];
	double max = time[0];

	for (int i = 0; i < times; i++) {
		//std::cout << i << "\ttime taken: " << time.at(i) << '\n';
		avg += time[i];
		if (min > time[i]) min = time[i];
		if (max < time[i]) max = time[i];
	}
	avg /= times;

	if (!isForCopy) {
		std::cout << "------TIME------"
			<< "\nAVG: " << avg
			<< "\nMIN: " << min
			<< "\nMAX: " << max << std::endl;
	}
	else {
		std::cout << avg << std::endl;
	}

	delete[] A;
	delete[] Y;
	delete[] X;
	delete[] time;
}

void testDecomposeSolveV2Par(int n, int times, bool isForCopy) {
	if (!isForCopy) {
		std::cout << "\nLU-decomposition (1ARR) parallel of matrix  " << n << " * " << n
			<< "\nRun " << times << " times." << std::endl;
	}

	double* A = new double[n * n];
	double* Y = new double[n];
	double* X = new double[n];
	double* time = new double[times];


	for (int i = 0; i < times; i++) {
		generateSystemOfLinearEquations(A, Y, n);
		/*init(L, n*n);
		init(U, n*n);
		init(Y, n*n);
		init(X, n*n);*/

		Timer t1;
		decompose_v2_par(A, n);
		solve(A, Y, X, n);
		time[i] = t1.elapsed();
	}

	double avg = 0.0;
	double min = time[0];
	double max = time[0];

	for (int i = 0; i < times; i++) {
		//std::cout << i << "\ttime taken: " << time.at(i) << '\n';
		avg += time[i];
		if (min > time[i]) min = time[i];
		if (max < time[i]) max = time[i];
	}
	avg /= times;

	if (!isForCopy) {
		std::cout << "------TIME------"
			<< "\nAVG: " << avg
			<< "\nMIN: " << min
			<< "\nMAX: " << max << std::endl;
	}
	else {
		std::cout << avg << std::endl;
	}

	delete[] A;
	delete[] Y;
	delete[] X;
	delete[] time;
}