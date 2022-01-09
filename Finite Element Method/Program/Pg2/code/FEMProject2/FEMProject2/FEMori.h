#pragma once

#include<iostream>
#include <iomanip>
#include <vector>
#include<cmath>
#include"ConjugateGradient.h"

using namespace std;

//d(x);
double d(double x){	return sin(x) + 2;}
//c(x);
double c(double x){	return x * x + 1;}
//真实解函数u(x)；
double u(double x){	return (x - 1) * sin(x);//return x * (x - 1);
}
//函数f(x)；
double f(double x){	return (-1 + x) * (1 + x * x) * sin(x) - cos(x) * ((-1 + x) * cos(x) + sin(x)) - (2 + sin(x)) * (2 * cos(x) - (-1 + x) * sin(x));
	//return -2 * (sin(x) + 2) - (2 * x - 1) * cos(x) + (x * x + 1) * x * (x - 1);
}

//第i个一次基函数；
double basis_1(int N, int i, double x)
{
	double he = 1.0 / N;
	if ((x > (i - 1.0) * he) && (x <= i * he))
		return (x - (i - 1.0) * he) / he;
	else if ((x > i * he) && (x < (i + 1.0) * he))
		return ((i + 1.0) * he - x) / he;
	else
		return 0.0;
}

//基于一次函数的复化三点Gauss积分公式，计算生成系数矩阵；
double gauss_integrate_k(int N, int i, int k, double left, double right)
{
	int n = 64;
	double delta_x = (right - left) / n;
	double temp = 0;
	if (i == k)
	{
		for (int j = 0; j < n; j++)
		{
			double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
			double x1 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 - sqrt(15) / 10 * delta_x;
			double x2 = (left + j * delta_x + left + (j + 1) * delta_x) / 2;
			double x3 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 + sqrt(15) / 10 * delta_x;
			temp += a1 * (N * N * d(x1) + c(x1) * basis_1(N, i, x1) * basis_1(N, i, x1)) + a2 * (N * N * d(x2) + c(x2) * basis_1(N, i, x2) * basis_1(N, i, x2)) + a3 * (N * N * d(x3) + c(x3) * basis_1(N, i, x3) * basis_1(N, i, x3));
		}
	}
	else
	{
		for (int j = 0; j < n; j++)
		{
			double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
			double x1 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 - sqrt(15) / 10 * delta_x;
			double x2 = (left + j * delta_x + left + (j + 1) * delta_x) / 2;
			double x3 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 + sqrt(15) / 10 * delta_x;
			temp += a1 * (-N * N * d(x1) + c(x1) * basis_1(N, i, x1) * basis_1(N, k, x1)) + a2 * (-N * N * d(x2) + c(x2) * basis_1(N, i, x2) * basis_1(N, k, x2)) + a3 * (-N * N * d(x3) + c(x3) * basis_1(N, i, x3) * basis_1(N, k, x3));
		}
	}
	return temp;
}

//基于一次函数基生成系数矩阵K;
void make_cof_k(int N, vector<vector<double>>& cof_k)
{
	//小区间个数有N个，系数矩阵cof_k大小为(N-1)*(N-1)；
	double he = 1.0 / N;//离散区间长度；
	for (int e = 1; e < N - 1; e++)
	{
		cof_k[e - 1][e - 1] = gauss_integrate_k(N, e, e, (e - 1) * he, (e + 1) * he);
		cof_k[e - 1][e] = gauss_integrate_k(N, e, e + 1, e * he, (e + 1) * he);
		cof_k[e][e - 1] = cof_k[e - 1][e];
	}
	int e = N - 1;
	cof_k[e - 1][e - 1] = gauss_integrate_k(N, e, e, (e - 1) * he, (e + 1) * he);
	/*
	cof_k[0][0] = (7.0 / 3 + he * he / 5 + (1 - cos(he)) / he) / he;
	for (int e = 2; e < N; e++)
	{
		double a = (e - 1.0) * he;
		cof_k[e - 2][e - 2] += ((7.0 / 3 + a * a / 3 + a * he / 6 + he * he / 30 + (cos(a) - cos(a + he)) / he) / he);
		cof_k[e - 1][e - 1] += ((7.0 / 3 + a * a / 3 + a * he / 2 + he * he / 5 + (cos(a) - cos(a + he)) / he) / he);
		cof_k[e - 2][e - 1] += ((-11.0 / 6 + a * a / 6 + a * he / 6 + he * he / 20 + (-cos(a) + cos(a + he)) / he) / he);
		cof_k[e - 1][e - 2] += ((-11.0 / 6 + a * a / 6 + a * he / 6 + he * he / 20 + (-cos(a) + cos(a + he)) / he) / he);
	}
	double a = (N - 1.0) * he;
	cof_k[N - 2][N - 2] += ((7.0 / 3 + a * a / 3 + a * he / 6 + he * he / 30 + (cos(a) - cos(a + he)) / he) / he);
	*/
}

//第i个二次基函数；
double basis_2(int N, int i, double x)
{
	double he = 1.0 / N;
	if (i < N)
	{
		if (x > (i - 1) * he && x <= i * he)
			return (2 * x - (2 * i - 1) * he) * (x - (i - 1) * he) / he / he;
		else if (x > i * he && x < (i + 1) * he)
			return (2 * x - (2 * i + 1) * he) * (x - (i + 1) * he) / he / he;
		else
			return 0;
	}
	else
	{
		if (x > (i - N) * he && x < (i - N + 1) * he)
			return 4 * (x - (i - N) * he) * ((i - N + 1) * he - x) / he / he;
		else
			return 0;
	}
}

//第i个二次基函数的导数；
double grad_basis_2(int N, int i, double x)
{
	double he = 1.0 / N;
	if (i < N)
	{
		if (x > (i - 1) * he && x <= i * he)
			return (4 * x - (4 * i - 3) * he) / he / he;
		else if (x > i * he && x < (i + 1) * he)
			return (4 * x - (4 * i + 3) * he) / he / he;
		else
			return 0;
	}
	else
	{
		if (x > (i - N) * he && x < (i - N + 1) * he)
			return 4 * ((2 * (i - N) + 1) * he - 2 * x) / he / he;
		else
			return 0;
	}
}

//基于二次函数的复化三点Gauss积分公式，计算生成系数矩阵；
double gauss_integrate_k_2(int N, int i, int k, double left, double right)
{
	int n = 2;
	double delta_x = (right - left) / n;
	double temp = 0;
	for (int j = 0; j < n; j++)
	{
		double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
		double x1 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 - sqrt(15) * delta_x / 10;
		double x2 = (left + j * delta_x + left + (j + 1) * delta_x) / 2;
		double x3 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 + sqrt(15) * delta_x / 10;
		temp += a1 * (grad_basis_2(N, i, x1) * grad_basis_2(N, k, x1) * d(x1) + c(x1) * basis_2(N, i, x1) * basis_2(N, k, x1))
			+ a2 * (grad_basis_2(N, i, x2) * grad_basis_2(N, k, x2) * d(x2) + c(x2) * basis_2(N, i, x2) * basis_2(N, k, x2))
			+ a3 * (grad_basis_2(N, i, x3) * grad_basis_2(N, k, x3) * d(x3) + c(x3) * basis_2(N, i, x3) * basis_2(N, k, x3));
	}
	return temp;
}

//基于二次函数基生成系数矩阵K；
void make_cof_k_2(int N, vector<vector<double>>& cof_k)
{
	//小区间个数有N个，系数矩阵cof_k大小为(2N-1)*(2N-1)；
	double he = 1.0 / N;//离散区间长度；
	/*
	for (int i = 1; i < N - 1; i++)
	{
		cof_k[i - 1][i - 1] = (4.0 / 15 + 28.0 / 3 / he / he + (2 * (24 * he + 8 * he * cos(he) + (-32 + he * he) * sin(he)) * sin(he * i)) / he / he / he / he / he) / he / he + (1.0 + 28.0 * i * i)/105;

		cof_k[i - 1][i] = (-1.0 / 30 + 2.0 / 3 / he / he + ((-32 + 3 * he * he) * cos(he * i) + (32 - 3 * he * he) * cos(he * (1 + i)) + 16 * he * (sin(he * i) + sin(he * (1 + i)))) / he / he / he / he / he) / he / he - (5.0 + 14.0 * i * (i + 1)) / 420;
		cof_k[i][i - 1] = cof_k[i - 1][i];

		cof_k[i - 1][i - 1 + N - 1] = (-4.0 * (-16.0 + he * he) * cos(i * he - he) / he / he / he / he / he + 4.0 * (-16.0 + 3 * he * he) * cos(i * he) / he / he / he / he / he - 16.0 / 3 / he / he + 1.0 / 15 - 24.0 * sin(he * i - he) / he / he / he / he - 40.0 * sin(he * i) / he / he / he / he) / he / he + (7.0 * i * i - 1.0) / 105;
		cof_k[i - 1 + N - 1][i - 1] = cof_k[i - 1][i - 1 + N - 1];

		cof_k[i - 1][i - 1 + N] = (-4.0 * (-16.0 + 3*he * he) * cos(i * he) / he / he / he / he / he + 4.0 * (-16.0 + he * he) * cos(i * he+he) / he / he / he / he / he - 16.0 / 3 / he / he + 1.0 / 15 - 24.0 * sin(he * i + he) / he / he / he / he - 40.0 * sin(he * i) / he / he / he / he) / he / he + (7.0 * i * i - 1.0) / 105;
		cof_k[i - 1 + N][i - 1] = cof_k[i - 1][i - 1 + N];

		cof_k[i - 1 + N - 1][i - 1 + N - 1] = (32.0 / 3 / he / he + 8.0 / 15 + 16.0 * (he * he - 8) * (cos(he * (i - 1)) - cos(he * (i - 1) + he)) / he / he / he / he / he + 64.0 * (sin(he * (i - 1)) + sin(he * (i - 1) + he)) / he / he / he / he) / he / he + (i - 1) * ((i - 1) - 1) * 1.0 / 105;
	}
	int i = N - 1;
	cof_k[N - 2][N - 2] = (4.0 / 15 + 28.0 / 3 / he / he + (2 * (24 * he + 8 * he * cos(he) + (-32 + he * he) * sin(he)) * sin(he * i)) / he / he / he / he / he) / he / he + (1.0 + 28.0 * i * i) / 105;

	cof_k[N - 1 - 1][N - 1 - 1 + N - 1] = (-4.0 * (-16.0 + he * he) * cos(i * he - he) / he / he / he / he / he + 4.0 * (-16.0 + 3 * he * he) * cos(i * he) / he / he / he / he / he - 16.0 / 3 / he / he + 1.0 / 15 - 24.0 * sin(he * i - he) / he / he / he / he - 40 * sin(he * i) / he / he / he / he) / he / he + (7.0 * i * i - 1.0) / 105;
	cof_k[N - 1 - 1 + N - 1][N - 1 - 1] = cof_k[N - 1 - 1][N - 1 - 1 + N - 1];

	cof_k[N - 1 - 1][N - 1 - 1 + N] = (-4.0 * (-16.0 + 3 * he * he) * cos(i * he) / he / he / he / he / he + 4.0 * (-16.0 + he * he) * cos(i * he + he) / he / he / he / he / he - 16.0 / 3 / he / he + 1.0 / 15 - 24.0 * sin(he * i + he) / he / he / he / he - 40.0 * sin(he * i) / he / he / he / he) / he / he + (7.0 * i * i - 1.0) / 105;
	cof_k[N - 1 - 1 + N][N - 1 - 1] = cof_k[N - 1 - 1][N - 1 - 1 + N];

	cof_k[2 * N - 3][2 * N - 3] = (32.0 / 3 / he / he + 8.0 / 15 + 16.0 * (he * he - 8) * (cos(he * (i - 1)) - cos(he * (i - 1) + he)) / he / he / he / he / he + 64.0 * (sin(he * (i - 1)) + sin(he * (i - 1) + he)) / he / he / he / he) / he / he + (i - 1) * ((i - 1) - 1) * 1.0 / 105;
	i = N;
	cof_k[2 * N - 2][2 * N - 2] = (32.0 / 3 / he / he + 8.0 / 15 + 16.0 * (he * he - 8) * (cos(he * (i - 1)) - cos(he * (i - 1) + he)) / he / he / he / he / he + 64.0 * (sin(he * (i - 1)) + sin(he * (i - 1) + he)) / he / he / he / he) / he / he + (i - 1) * ((i - 1) - 1) * 1.0 / 105;
	*/

	
	for (int i = 1; i < N - 1; i++)
	{
		cof_k[i - 1][i - 1] = gauss_integrate_k_2(N, i, i, (i - 1) * he, (i + 1) * he);

		cof_k[i - 1][i] = gauss_integrate_k_2(N, i, i + 1, i * he, (i + 1) * he);
		cof_k[i][i - 1] = cof_k[i - 1][i];

		cof_k[i - 1][i - 1 + N - 1] = gauss_integrate_k_2(N, i, i + N - 1, (i - 1) * he, i * he);
		cof_k[i - 1 + N - 1][i - 1] = cof_k[i - 1][i - 1 + N - 1];

		cof_k[i - 1][i - 1 + N] = gauss_integrate_k_2(N, i, i + N, i * he, (i + 1) * he);
		cof_k[i - 1 + N][i - 1] = cof_k[i - 1][i - 1 + N];

		cof_k[i - 1 + N - 1][i - 1 + N - 1] = gauss_integrate_k_2(N, i + N - 1, i + N - 1, (i - 1) * he, i * he);
	}
	int i = N - 1;
	cof_k[N - 2][N - 2] = gauss_integrate_k_2(N, i, i, (i - 1) * he, (i + 1) * he);

	cof_k[N - 1 - 1][N - 1 - 1 + N - 1] = gauss_integrate_k_2(N, i, i + N - 1, (i - 1) * he, i * he);
	cof_k[N - 1 - 1 + N - 1][N - 1 - 1] = cof_k[N - 1 - 1][N - 1 - 1 + N - 1];

	cof_k[N - 1 - 1][N - 1 - 1 + N] = gauss_integrate_k_2(N, i, i + N, i * he, (i + 1) * he);
	cof_k[N - 1 - 1 + N][N - 1 - 1] = cof_k[N - 1 - 1][N - 1 - 1 + N];

	cof_k[2 * N - 3][2 * N - 3] = gauss_integrate_k_2(N, i + N - 1, i + N - 1, (i - 1) * he, i * he);

	i = N;
	cof_k[2 * N - 2][2 * N - 2] = gauss_integrate_k_2(N, i + N - 1, i + N - 1, (i - 1) * he, i * he);
	
	/*
	cof_k[0][0] = a[1][1] / he;
	cof_k[0][1] = a[1][2] / he;
	cof_k[1][0] = a[2][1] / he;
	cof_k[1][1] = a[2][2] / he;
	int e = 3;
	while (e < 2 * N - 2)
	{
		cof_k[index_2(e, 0)][index_2(e, 0)] += (a[0][0] / he);
		cof_k[index_2(e, 0)][index_2(e, 1)] += (a[0][1] / he);
		cof_k[index_2(e, 0)][index_2(e, 2)] += (a[0][2] / he);
		cof_k[index_2(e, 1)][index_2(e, 0)] += (a[1][0] / he);
		cof_k[index_2(e, 1)][index_2(e, 1)] += (a[1][1] / he);
		cof_k[index_2(e, 1)][index_2(e, 2)] += (a[1][2] / he);
		cof_k[index_2(e, 2)][index_2(e, 0)] += (a[2][0] / he);
		cof_k[index_2(e, 2)][index_2(e, 1)] += (a[2][1] / he);
		cof_k[index_2(e, 2)][index_2(e, 2)] += (a[2][2] / he);
		e += 2;
	}
	/*for (int e = 2; e < N; e++)
	{
		cof_k[index(e, 0)][index(e, 0)] += (a[0][0] / he);
		cof_k[index(e, 0)][index(e, 1)] += (a[0][1] / he);
		cof_k[index(e, 0)][index(e, 2)] += (a[0][2] / he);
		cof_k[index(e, 1)][index(e, 0)] += (a[1][0] / he);
		cof_k[index(e, 1)][index(e, 1)] += (a[1][1] / he);
		cof_k[index(e, 1)][index(e, 2)] += (a[1][2] / he);
		cof_k[index(e, 2)][index(e, 0)] += (a[2][0] / he);
		cof_k[index(e, 2)][index(e, 1)] += (a[2][1] / he);
		cof_k[index(e, 2)][index(e, 2)] += (a[2][2] / he);
	}
	cof_k[2 * N - 3][2 * N - 3] += a[0][0] / he;
	cof_k[2 * N - 3][2 * N - 2] += a[0][1] / he;
	cof_k[2 * N - 2][2 * N - 3] += a[1][0] / he;
	cof_k[2 * N - 2][2 * N - 2] += a[1][1] / he;
	*/
}

//基于一次函数的复化三点Gauss积分公式；
double gauss_integrate_1(int N, int i, double left, double right)
{
	int n = 2;
	double delta_x = (right - left) / n;
	double temp = 0;
	for (int j = 0; j < n; j++)
	{
		double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
		double x1 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 - sqrt(15) / 10 * delta_x;
		double x2 = (left + j * delta_x + left + (j + 1) * delta_x) / 2;
		double x3 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 + sqrt(15) / 10 * delta_x;
		temp += a1 * f(x1) * basis_1(N, i, x1) + a2 * f(x2) * basis_1(N, i, x2) + a3 * f(x3) * basis_1(N, i, x3);
	}
	return temp;
}

//基于二次函数的复化三点Gauss积分公式；
double gauss_integrate_2(int N, int i, double left, double right)
{
	int n = 128;
	double delta_x = (right - left) / n;
	double temp = 0;
	for (int j = 0; j < n; j++)
	{
		double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
		double x1 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 - sqrt(15) / 10 * delta_x;
		double x2 = (left + j * delta_x + left + (j + 1) * delta_x) / 2;
		double x3 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 + sqrt(15) / 10 * delta_x;
		temp += a1 * f(x1) * basis_2(N, i, x1) + a2 * f(x2) * basis_2(N, i, x2) + a3 * f(x3) * basis_2(N, i, x3);
	}
	return temp;
}

//生成右侧向量；
void make_cof_f(int N, vector<double>& cof_f_1, vector<double>& cof_f_2)
{
	double he = 1.0 / N;
	for (int i = 1; i < N; i++)
	{
		cof_f_1[i - 1] = gauss_integrate_1(N, i, (i - 1) * he, (i + 1) * he);
		cof_f_2[i - 1] = gauss_integrate_2(N, i, (i - 1) * he, (i + 1) * he);
	}
	for (int i = N; i < 2 * N; i++)
	{
		cof_f_2[i - 1] = gauss_integrate_2(N, i, (i - N) * he, (i - N + 1) * he);
	}
}

//基于1次基函数的求解函数uh(x)；
double uh_basis_1(int N, vector<double> cof_u, double x)
{
	double temp = 0;
	double he = 1.0 / N;
	for (int i = 1; i < N; i++)
		temp += cof_u[i - 1] * basis_1(N, i, x);
	return temp;
}

//离散计算无穷范数误差；
double error_max_1(int N, vector<double> cof_u)
{
	int n = 100;
	double delta_x = 1.0 / n;
	vector<double> error(n + 1, 0);
	//cout << endl << "N=" << N << "数值解" << endl;
	for (int i = 0; i <= n; i++)
		//cout << "{" << i * delta_x << "," << uh_basis_1(N, cof_u, i * delta_x) << "},",
		error[i] = fabs(u(i * delta_x) - uh_basis_1(N, cof_u, i * delta_x));
	double max = error[0];
	for (int i = 1; i <= n; i++)
		if (error[i] > max)
			max = error[i];
	return max;
}

//基于复化Gauss积分公式离散计算L1范数误差；
double error_L1_1(int N, vector<double> cof_u)
{
	int n = 128;
	double delta_x = 1.0 / n;
	double temp = 0;
	for (int j = 0; j < n; j++)
	{
		double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
		double x1 = (j * delta_x + (j + 1) * delta_x) / 2 - sqrt(15) / 10 * delta_x;
		double x2 = (j * delta_x + (j + 1) * delta_x) / 2;
		double x3 = (j * delta_x + (j + 1) * delta_x) / 2 + sqrt(15) / 10 * delta_x;
		temp += a1 * fabs(u(x1) - uh_basis_1(N, cof_u, x1)) + a2 * fabs(u(x2) - uh_basis_1(N, cof_u, x2)) + a3 * fabs(u(x3) - uh_basis_1(N, cof_u, x3));
	}
	return temp;
}

//基于2次基函数的求解函数uh(x)；
double uh_basis_2(int N, vector<double> cof_u, double x)
{
	double temp = 0.0;
	double he = 1.0 / N;
	for (int i = 1; i < 2 * N; i++)
		temp += cof_u[i - 1] * basis_2(N, i, x);
	return temp;
}

//离散计算无穷范数误差；
double error_max_2(int N, vector<double> cof_u)
{
	int n = 4*N;
	double delta_x = 1.0 / n;
	vector<double> error(n + 1, 0);
	for (int i = 0; i <= n; i++)
		error[i] = fabs(u(i * delta_x) - uh_basis_2(N, cof_u, i * delta_x));
	double max = error[0];
	for (int i = 1; i <= n; i++)
		if (error[i] > max)
			max = error[i];
	return max;
}

//基于复化Gauss积分公式离散计算L1范数误差；
double error_L1_2(int N, vector<double> cof_u)
{
	int n = 128;
	double delta_x = 1.0 / n;
	double temp = 0;
	for (int j = 0; j < n; j++)
	{
		double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
		double x1 = (j * delta_x + (j + 1) * delta_x) / 2 - sqrt(15) / 10 * delta_x;
		double x2 = (j * delta_x + (j + 1) * delta_x) / 2;
		double x3 = (j * delta_x + (j + 1) * delta_x) / 2 + sqrt(15) / 10 * delta_x;
		temp += a1 * fabs(u(x1) - uh_basis_2(N, cof_u, x1)) + a2 * fabs(u(x2) - uh_basis_2(N, cof_u, x2)) + a3 * fabs(u(x3) - uh_basis_2(N, cof_u, x3));
	}
	return temp;
}
void vector_output(vector<double> vec)
{
	//cout << "向量输出如下：" << endl;
	for (int i = 0; i < vec.size(); i++)
		cout << vec[i] << "\t";
	cout << endl;
}
void matrix_output(vector<vector<double>> matrix)
{
	//cout << endl << "矩阵输出如下" << endl;
	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix[0].size(); j++)
			cout << matrix[i][j] << "\t";
		cout << endl;
	}
	cout << endl;
}
//基函数计算求解，输出无穷范数误差；
void error_out(vector<vector<double>>& error)
{
	for (int i = 0; i < 5; i++)
	{
		//区间离散程度；
		int N = pow(2, i) * 10;
		double he = 1.0 / N;
		//系数矩阵K；
		vector<vector<double>> cof_k_1(N - 1, vector<double>(N - 1, 0)), cof_k_2(2 * N - 1, vector<double>(2 * N - 1, 0));
		make_cof_k(N, cof_k_1);
		make_cof_k_2(N, cof_k_2);
		//右侧函数f与基内积向量；
		vector<double> cof_f_1(N - 1, 0), cof_f_2(2 * N - 1, 0);
		make_cof_f(N, cof_f_1, cof_f_2);
		//离散解的系数向量；
		vector<double> cof_u_1(N - 1, 0), cof_u_2(2 * N - 1, 1);
		//共轭梯度法求解；
		int times = 0;             //共轭梯度迭代法迭代次数；
		double iteration_error = 0;//共轭梯度迭代法的输出参考误差；
		Conjugate_Gradient(N - 1, cof_k_1, cof_f_1, 1e-20, cof_u_1, times, iteration_error);
		times = 0;           //共轭梯度迭代法迭代次数；
		iteration_error = 0;//共轭梯度迭代法的输出参考误差；
		Conjugate_Gradient(2 * N - 1, cof_k_2, cof_f_2, 1e-20, cof_u_2, times, iteration_error);
		if (i == 0)
			matrix_output(cof_k_2),
			cout << endl,
			vector_output(cof_f_2),
			cout << endl,
			vector_output(cof_u_2);
		//计算无穷范数作为误差；
		error[i][0] = error_max_1(N, cof_u_1);
		error[i][2] = error_L1_1(N, cof_u_1);
		error[i][4] = error_max_2(N, cof_u_2);
		error[i][6] = error_L1_2(N, cof_u_2);
		if (i > 0)
		{
			error[i][1] = log(error[i - 1][0] / error[i][0]) / log(2);
			error[i][3] = log(error[i - 1][2] / error[i][2]) / log(2);
			error[i][5] = log(error[i - 1][4] / error[i][4]) / log(2);
			error[i][7] = log(error[i - 1][6] / error[i][6]) / log(2);
		}
	}
}