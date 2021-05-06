#pragma once

void swap(double& a, double& b);
void swap(int& a, int& b);
void copy(double* from, double* to, int rows, int cols);
void copy(double* from, double* to, int size);
void print(double* a, double* y, int n);
void print(double* a, int size);
double maxInMatrix(double* matrix, int size);
int getCountsOfDigits(int number);
void printMatrix(double* matrix, int size);
void generateSystemOfLinearEquations(double* a, double* y, int n);
void generateRandomMatrix(double* matrix, int size);
void init(double* a, int size, double elem);
void init(double* a, int size);
void input(double* a, double* y, int n);