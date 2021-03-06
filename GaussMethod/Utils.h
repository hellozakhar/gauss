#pragma once
#include <iostream>

void swap(double& a, double& b) {
	double tmp = a;
	a = b;
	b = tmp;
}

void swap(int& a, int& b) {
	int tmp = a;
	a = b;
	b = tmp;
}

void copy(double* from, double* to, int rows, int cols) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			to[i * rows + j] = from[i * rows + j];
		}
	}
}

void copy(double* from, double* to, int size) {
	for (int i = 0; i < size; i++) {
		to[i] = from[i];
	}
}

void print(double* a, double* y, int n) {
	std::cout << std::endl;
	char sep = ' ';
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			//std::cout << a[i * n + j] << "*x" << j;
			std::cout << a[i * n + j];
			if (j < n - 1)
				std::cout << sep;
		}
		std::cout << " = " << y[i] << std::endl;
	}
}

// TO-DO: make check "determinant == 0 ?"
// TO-DO: change the method of generating numbers
void generateLinearSystem(double* a, double* y, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			a[i * n + j] = rand() % 9999;
		}

		y[i] = rand() % 9999;
	}
}

void init(double* a, int size, double elem) {
	for (int i = 0; i < size; i++) {
		a[i] = elem;
	}
}

void init(double* a, int size) {
	init(a, size, 0.0);
}