#pragma once

void gaussDef(double* a, double* y, double* x, int n);
void gaussMax(double* a, double* y, double* x, int n);
void decompose_v1(double* A, double* L, double* U, int n);
void solve(double* L, double* U, double* y, double* x, int n);
void solve(double* LU, double* y, double* x, int n);
void decompose_v2(double* A, int n);
void decompose_v2_par(double* A, int n);