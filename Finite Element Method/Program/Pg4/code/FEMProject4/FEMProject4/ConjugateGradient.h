#pragma once

#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

//矩阵乘以列向量；
void Matrix_times_vector(vector<vector<double>> matrix, vector<double> column, vector<double>& value)
{
	int dim = size(matrix);
	for (int i = 0; i < dim; i++)
	{
		value[i] = 0;
		for (int j = 0; j < dim; j++)
			value[i] += matrix[i][j] * column[j];
	}
}

//向量之间作减法；
void Vector_subtraction(vector<double> x, vector<double> y, vector<double>& dx)
{
	int dim = 0;
	dim = size(x);
	for (int i = 0; i < dim; i++)
		dx[i] = x[i] - y[i];
}

//向量作内积；
double Inner_product_of_vectors(vector<double> x, vector<double> y)
{
	int i, n;
	double temp = 0;
	n = size(x);
	for (i = 0; i < n; i++)
		temp += x[i] * y[i];
	return temp;
}

//向量的2-范数，并返回范数值；
double vector_2_norm(vector<double> x)
{
	int dim = size(x);
	double temp = 0.0;
	for (int i = 0; i < dim; i++)
		temp = temp + x[i] * x[i];
	return sqrt(temp);
}

//辅助函数ptAp；
double pAp(vector<vector<double>> matrix, vector<double> column)
{
	int dim = size(matrix);
	vector<double> value(dim, 0);
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			value[i] += matrix[i][j] * column[j];
	return Inner_product_of_vectors(column, value);
}

//辅助函数，共轭梯度迭代生成新的xk；
void xk(vector<double>& x, int dim, vector<double> x0, double a0, vector<double> p0)
{
	for (int i = 0; i < dim; i++)
		x[i] = x0[i] + a0 * p0[i];
}

//辅助函数，生成残向量rk；
void rk(vector<double>& r, int dim, vector<double> r0, vector<double> p0, double a0, vector<vector<double>> matrix)
{
	vector<double> temp(dim, 0);
	Matrix_times_vector(matrix, p0, temp);
	for (int i = 0; i < dim; i++)
		r[i] = r0[i] - a0 * temp[i];
}

//共轭梯度法；
void Conjugate_Gradient(int dim, vector<vector<double>> matrix, vector<double> b, double epsilon, vector<double>& x, int& times, double& error)
{
	double e = epsilon * vector_2_norm(b);
	vector<double> x0(dim, 0);
	//x0[4] = 1;
	vector<double> r0(dim, 0), r(dim, 0), p0(dim, 0);
	Matrix_times_vector(matrix, x0, r0);
	Vector_subtraction(b, r0, r0);//顺序?
	r = r0;
	x = x0;
	times = 0;
	while (vector_2_norm(r) > e && times < 10000)
	{
		times++;
		if (times == 1) p0 = r0;
		else
		{
			double beta = Inner_product_of_vectors(r, r) / Inner_product_of_vectors(r0, r0);
			xk(p0, dim, r, beta, p0);
		}
		double alpha = Inner_product_of_vectors(r, r) / pAp(matrix, p0);
		xk(x, dim, x, alpha, p0);
		r0 = r;//储存rk-1;
		rk(r, dim, r0, p0, alpha, matrix);
	}
	error = vector_2_norm(r);
}
