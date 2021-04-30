#pragma once
#include <iostream>
#include <iomanip>

#include "GaussMethodAux.h"

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

// Prints linear equation system to console
void print(double* a, double* y, int n) {
	std::cout << std::endl;
	char sep = ' ';
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << a[i * n + j];
			if (j < n - 1)
				std::cout << sep;
		}
		std::cout << " = " << y[i] << std::endl;
	}
}

// Prints content of array to console
void print(double* a, int size) {
	char sep = ' ';
	for (int i = 0; i < size; i++) {
		std::cout << a[i] << sep;
	}
	std::cout << std::endl;
}

/*
 * Функция maxInMatrix возвращает наибольший элемент матрицы matrix размера size.
 * см. проект MultMatrix.
 */
double maxInMatrix(double* matrix, int size) {
	double max = matrix[0];
	for (int i = 0; i < size*size; i++) {
		if (max < matrix[i]) {
			max = matrix[i];
		}
	}
	return max;
}

/*
 * Функция getCountsOfDigits возвращает количество цифр целого числа number.
 * см. проект MultMatrix.
 */
int getCountsOfDigits(int number) {
	int count = (number == 0) ? 1 : 0;
	while (number != 0) {
		count++;
		number /= 10;
	}
	return count;
}

/*
 * Функция printMatrix выводит матрицу matrix размера size в консоль.
  * см. проект MultMatrix.
 */
void printMatrix(double* matrix, int size) {
	int width_elem = (int)getCountsOfDigits(maxInMatrix(matrix, size)) + 1;
	std::cout << '\n';
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			std::cout.width(width_elem);
			std::cout << matrix[i*size + j] << " ";
		}
		std::cout << '\n';
	}
}

// TO-DO: make check "determinant == 0 ?"
// TO-DO: change the method of generating numbers
void generateSystemOfLinearEquations(double* a, double* y, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			a[i * n + j] = rand() % 9999;
		}

		y[i] = rand() % 9999;
	}
}

/*
* Функция generateRandomMatrix заполняет матрицу matrix
* размера size случайно сгенерированными элементами типа double.
* см. проект MultMatrix.
*/
void generateRandomMatrix(double* matrix, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			matrix[i*size + j] = rand() % 9999;
		}
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

void input(double* a, double* y, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << "a[" << i << "][" << j << "]= ";
			std::cin >> a[i * n + j];
		}
	}

	for (int i = 0; i < n; i++) {
		std::cout << "y[" << i << "]= ";
		std::cin >> y[i];
	}
}