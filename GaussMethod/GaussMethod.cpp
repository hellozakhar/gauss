#include <iostream>
#include <omp.h>

#include "GaussMethodAux.h"
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

void decompose_v1(double* A, double* L, double* U, int n) {
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

void solve(double* L, double* U, double* y, double* x, int n) {
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

void solve(double* LU, double* y, double* x, int n) {
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

void decompose_v2(double* A, int n) {
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

void decompose_v2_par(double* A, int n) {
	for (int i = 0; i < n; i++) {

		#ifdef _OPENMP
		omp_set_num_threads(4);
		#pragma omp parallel for
		#endif
		for (int k = i + 1; k < n; k++) {
			double mu = A[k * n + i] / A[i * n + i];

			for (int j = i; j < n; j++) {
				A[k * n + j] -= mu * A[i * n + j];
			}

			A[k * n + i] = mu;
		}
	}
}