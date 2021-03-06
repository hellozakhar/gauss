#include <iostream>

#include "Utils.h"
#include "Timer.h"
#include "GaussMethod.h"

const double EPS = 1e-5;

void gaussDef(double* a, double* y, double* x, int n) {

	/*  Прямой ход метода Гаусса.
		На итерации i, 1 <= i <= n, метода производится исключение
		неизвестной i для всех уравнений с номерами k, i < k <= n.
		Для этого из этих уравнений осуществляется вычитание строки i,
		умноженной на константу(aki / aii), чтобы результирующий
		коэффициент при неизвестной xi в строках оказался
		нулевым.
		*/

	// Множитель Гаусса
	double mu;
	
	for (int i = 0; i < n - 1; i++) {
		for (int k = i + 1; k <= n - 1; k++) {
			
			mu = a[k * n + i] / a[i * n + i];

			for (int j = i; j <= n - 1; j++) {
				a[k * n +	j] = a[k * n + j] - mu * a[i * n + j];
				
			}

			y[k] = y[k] - mu * y[i];
		}
	}

	/*  Обратный ход метода Гаусса.
		После приведения матрицы коэффициентов к треугольному
		виду становится возможным определение значений неизвестных :
		- Из последнего уравнения преобразованной системы может
		быть вычислено значение переменной xn,
		- Из предпоследнего уравнения становится возможным
		определение переменной xn - 1, и т.д.
		*/

	for (int i = n - 1; i >= 0; i--) {
		x[i] = y[i];
		
		for (int j = i + 1; j < n; j++) {
			x[i] -= a[i * n + j] * x[j];
		}

		x[i] /= a[i * n + i];
	}
}

void gaussMax(double* a, double* y, double* x, int n) {
	double max;
	double mu; // Множитель Гаусса
	int index;
	int* row = new int[n];

	for (int i = 0; i < n; i++) {
		row[i] = i;
	}

	/*  Прямой ход метода Гаусса.
		На итерации i, 1 <= i <= n, метода производится исключение
		неизвестной i для всех уравнений с номерами k, i < k <= n.
		Для этого из этих уравнений осуществляется вычитание строки i,
		умноженной на константу(aki / aii), чтобы результирующий
		коэффициент при неизвестной xi в строках оказался
		нулевым.
		*/

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

	/*  Обратный ход метода Гаусса.
		После приведения матрицы коэффициентов к треугольному
		виду становится возможным определение значений неизвестных :
		- Из последнего уравнения преобразованной системы может
		быть вычислено значение переменной xn,
		- Из предпоследнего уравнения становится возможным
		определение переменной xn - 1, и т.д.
		*/

	for (int i = n - 1; i >= 0; i--) {
		x[i] = y[row[i]];

		for (int j = i + 1; j < n; j++) {
			x[i] -= a[row[i] * n + j] * x[j];
		}

		x[i] /= a[row[i] * n + i];
	}	 
}

void testGaussDef(bool isAuto) {
	int n;
	std::cout << "Enter count of equations: ";
	std::cin >> n;

	double* a = new double[n * n];
	double* y = new double[n];
	double* x = new double[n];
	init(x, n);

	if (isAuto) {
		generateLinearSystem(a, y, n);
	}
	else {
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
		generateLinearSystem(a, y, n);
	}
	else {
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

	gaussMax(a, y, x, n);

	for (int i = 0; i < n; i++) {
		std::cout << "x[" << i << "]=" << x[i] << std::endl;
	}

	delete[] a;
	delete[] y;
	delete[] x;
}

int main() {
	// isAuto true  => coeffs and y - generated,
	// isAuto false => coeffs and y - enter manually.
	const bool isAuto = true;

	testGaussDef(isAuto);

	std::cin.get();
	return 0;
}