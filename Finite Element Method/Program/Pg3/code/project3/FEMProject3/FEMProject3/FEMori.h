#pragma once

#include<iostream>
#include <iomanip>
#include <vector>
#include<cmath>

#include"ConjugateGradient.h"

using namespace std;
#define epsilon7 1e-7
//d(x);
double d(double x)
{
	double epsilon = 1e-1;
	return epsilon;
}

//c(x);
double c(double x)
{
	return 1.0;
}

//真实解函数u(x)；
double u(double x)
{
	return 0.6 * (1.0 - 1.0 * exp(10.0 * x)) / (exp(10.0) - 1.0) + 0.5 * x * x + 0.1 * x;
}

//epsilon=1e-7 真实解函数u(x)；
double u_2(double x)
{
	double c1 = epsilon7 * x;
	double c0 = x * x / 2;
	double c2 = (x - 1.0) / epsilon7;
	return c0 + c1 - (epsilon7 + 0.5) * exp(c2);
}

//函数f(x)；
double f(double x)
{
	return x;
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

//第i个基函数的导数；
double basis_diff_1(int N, int i, double x)
{
	double he = 1.0 / N;
	if ((x > (i - 1.0) * he) && (x <= i * he))
		return 1.0 / he;
	else if ((x > i * he) && (x < (i + 1.0) * he))
		return -1.0 / he;
	else
		return 0.0;
}

//基于一次函数的复化三点Gauss积分公式，计算生成系数矩阵；
double gauss_integrate_k(int N, int i, int k, double left, double right)
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
		temp += a1 * (basis_diff_1(N, i, x1) * basis_diff_1(N, k, x1) * d(x1) + c(x1) * basis_1(N, i, x1) * basis_diff_1(N, k, x1))
			+ a2 * (basis_diff_1(N, i, x2) * basis_diff_1(N, k, x2) * d(x2) + c(x2) * basis_1(N, i, x2) * basis_diff_1(N, k, x2))
			+ a3 * (basis_diff_1(N, i, x3) * basis_diff_1(N, k, x3) * d(x3) + c(x3) * basis_1(N, i, x3) * basis_diff_1(N, k, x3));
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
		cof_k[e][e - 1] = gauss_integrate_k(N, e + 1, e, e * he, (e + 1) * he);
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

//生成一次空间的右侧向量；
void make_cof_f(int N, vector<double>& cof_f_1)
{
	double he = 1.0 / N;
	for (int i = 1; i < N; i++)
		cof_f_1[i - 1] = gauss_integrate_1(N, i, (i - 1) * he, (i + 1) * he);
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
	int n = 1000;
	double delta_x = 1.0 / n;
	vector<double> error(n + 1, 0);
	//cout << endl << "N=" << N << "数值解" << endl;
	for (int i = 0; i <= n; i++)
	{
		/*
		int k = 0;
		double temp = fabs(u(i * delta_x) - uh_basis_1(N, cof_u, i * delta_x));
		if (temp!=0)
			while (temp <= 1)
			{
				k++;
				temp *= 10;
			}
		cout << "{" << i * delta_x << "," << temp << "*10^(" << -k << ")" << "},";
		*/
		error[i] = fabs(u(i * delta_x) - uh_basis_1(N, cof_u, i * delta_x));
	}
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

//三对角阵Gauss消去求解方程组；
void gauss(vector<vector<double>>& cof_k, vector<double>& cof_f, vector<double>& cof_u)
{
	//消去下三角并修改右侧系数；
	for (int i = 1; i < cof_k.size(); i++)
	{
		cof_k[i][i] = cof_k[i][i] - cof_k[i][i - 1] / cof_k[i - 1][i - 1] * cof_k[i - 1][i];
		cof_f[i] = cof_f[i] - cof_k[i][i - 1] / cof_k[i - 1][i - 1] * cof_f[i - 1];
		cof_k[i][i - 1] = 0;
	}

	//求解方程组；
	for (int i = cof_k.size() - 1; i > 0; i--)
	{
		cof_u[i] = cof_f[i] / cof_k[i][i];
		//cout << cof_f[i] <<"  "<< cof_k[i][i] <<"  "<< cof_u[i] << endl;
		cof_f[i - 1] -= cof_u[i] * cof_k[i - 1][i];
	}
	cof_u[0] = cof_f[0] / cof_k[0][0];
}

/*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Shinshkin下构造一次函数基；
double basis_shinshkin_1(int N, int i, double left, double right, double x)
{
	double he = (right - left) / N;
	if ((x > left + (i - 1.0) * he) && (x <= left + i * he))
		return (x - left - (i - 1.0) * he) / he;
	else if ((x > left + i * he) && (x < left + (i + 1.0) * he))
		return (left + (i + 1.0) * he - x) / he;
	else
		return 0.0;
}

//Shinshkin下构造一次函数基的导数；
double basis_shinshkin_diff_1(int N, int i, double left, double right, double x)
{
	double he = (right - left) / N;
	if ((x > left + (i - 1.0) * he) && (x <= left + i * he))
		return 1.0 / he;
	else if ((x > left + i * he) && (x < left + (i + 1.0) * he))
		return -1.0 / he;
	else
		return 0.0;
}

//Shinshkin下构造系数矩阵
double gauss_integrate_shinshkin_k(int N, double basisleft, double basisright, int i, int k, double left, double right)
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
		temp += a1 * (basis_shinshkin_diff_1(N, i, basisleft, basisright, x1) * basis_shinshkin_diff_1(N, k, basisleft, basisright, x1) * d(x1) + c(x1) * basis_shinshkin_1(N, i, basisleft, basisright, x1) * basis_shinshkin_diff_1(N, k, basisleft, basisright, x1))
			+ a2 * (basis_shinshkin_diff_1(N, i, basisleft, basisright, x2) * basis_shinshkin_diff_1(N, k, basisleft, basisright, x2) * d(x2) + c(x2) * basis_shinshkin_1(N, i, basisleft, basisright, x2) * basis_shinshkin_diff_1(N, k, basisleft, basisright, x2))
			+ a3 * (basis_shinshkin_diff_1(N, i, basisleft, basisright, x3) * basis_shinshkin_diff_1(N, k, basisleft, basisright, x3) * d(x3) + c(x3) * basis_shinshkin_1(N, i, basisleft, basisright, x3) * basis_shinshkin_diff_1(N, k, basisleft, basisright, x3));
	}
	return temp;
}

//Shinshkin下基于一次函数基生成系数矩阵K;
void make_cof_shinshkin_k(int N, double epsilon, vector<vector<double>>& cof_k)
{
	double tao = -2 * epsilon * log(N) + 1;
	//小区间个数有N个，系数矩阵cof_k大小为(N-1)*(N-1)；
	double he_1 = tao / N;//离散区间长度；
	double he_2 = (1.0 - tao) / N;
	
	/*
	for (int e = 1; e < N - 1; e++)
	{
		cof_k_1[e - 1][e - 1] = gauss_integrate_shinshkin_k(N, 0, tao, e, e, (e - 1) * he_1, (e + 1) * he_1);
		cof_k_1[e - 1][e] = gauss_integrate_shinshkin_k(N, 0, tao, e, e + 1, e * he_1, (e + 1) * he_1);
		cof_k_1[e][e - 1] = gauss_integrate_shinshkin_k(N, 0, tao, e + 1, e, e * he_1, (e + 1) * he_1);

		cof_k_2[e - 1][e - 1] = gauss_integrate_shinshkin_k(N, tao, 1.0-tao, e, e, tao+(e - 1) * he_2, tao+(e + 1) * he_2);
		cof_k_2[e - 1][e] = gauss_integrate_shinshkin_k(N, tao, 1.0 - tao, e, e + 1, tao+e * he_2, tao+(e + 1) * he_2);
		cof_k_2[e][e - 1] = gauss_integrate_shinshkin_k(N, tao, 1.0 - tao, e + 1, e, tao+e * he_2, tao+(e + 1) * he_2);
	}
	int e = N - 1;
	cof_k_1[e - 1][e - 1] = gauss_integrate_shinshkin_k(N, 0, tao, e, e, (e - 1) * he_1, (e + 1) * he_1);
	cof_k_2[e - 1][e - 1] = gauss_integrate_shinshkin_k(N, tao, 1.0 - tao, e, e, tao+(e - 1) * he_2, tao+(e + 1) * he_2);
	//
	///*
	for (int e = 1; e < N - 1; e++)
	{
		cof_k[e - 1][e - 1] = epsilon * 2 / he_1;
		cof_k[e - 1][e] = -epsilon / he_1 + 0.5;
		cof_k[e][e - 1] = -epsilon / he_1 - 0.5;

		cof_k[e - 1 + N][e - 1 + N] = epsilon * 2 / he_2;
		cof_k[e - 1 + N][e - 1 + N + 1] = -epsilon / he_2 + 0.5;
		cof_k[e - 1 + N + 1][e - 1 + N] = -epsilon / he_2 - 0.5;
	}
	int e = N - 1;
	cof_k[e - 1][e - 1] = epsilon * 2 / he_1;
	cof_k[e - 1][e] = -epsilon / he_1 + 0.5;
	cof_k[e][e - 1] = -epsilon / he_1 - 0.5;

	cof_k[e - 1 + N][e - 1 + N] = epsilon * 2 / he_2;

	e = N;
	cof_k[e - 1][e - 1] = epsilon / he_1 + epsilon / he_2;
	cof_k[e - 1][e] = -epsilon / he_2 + 0.5;
	cof_k[e][e - 1] = -epsilon / he_2 - 0.5;

	//
}

//Shinshkin下基于一次函数的复化三点Gauss积分公式；
double gauss_integrate_shinshkin_1(int N, double basisleft, double basisright, int i, double left, double right)
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
		temp += a1 * f(x1) * basis_shinshkin_1(N, i, basisleft, basisright, x1) + a2 * f(x2) * basis_shinshkin_1(N, i, basisleft, basisright, x2) + a3 * f(x3) * basis_shinshkin_1(N, i, basisleft, basisright, x3);
	}
	return temp;
}

//生成一次空间的右侧向量；
void make_cof_shinshkin_f(int N, double tao, vector<double>& cof_f)
{
	double he_1 = tao / N;
	double he_2 = (1.0 - tao) / N;
	for (int i = 1; i < N; i++)
	{
		cof_f[i - 1] = i * he_1 * he_1;
		cof_f[i - 1 + N] = (i * he_2 + tao) * he_2;// + he_2 * tao;
		//cof_f_1[i - 1] = gauss_integrate_shinshkin_1(N, 0, tao, i, (i - 1) * he_1, (i + 1) * he_1);
		//cof_f_2[i - 1] = gauss_integrate_shinshkin_1(N, tao, 1 - tao, i, tao+(i - 1) * he_2, tao+(i + 1) * he_2);
	}
	int i = N;
	cof_f[i - 1] = he_1 / 2 + (-2 * tao + 1.0) / 6 / N / N;
}

//基于1次基函数的求解函数uh(x)；
double uh_basis_shinshkin_1(int N, double tao, vector<double> cof_u, double x)
{
	double he_1 = tao / N;
	double he_2 = (-tao + 1.0) / N;
	double temp = 0;
	for (int i = 1; i < N; i++)
		//cout<<endl<< cof_u_1[i - 1] * basis_shinshkin_1(N, i, 0, tao, x) + cof_u_2[i - 1] * basis_shinshkin_1(N, i, tao, 1 - tao, x),
		temp += (cof_u[i - 1] * basis_shinshkin_1(N, i, 0, tao, x) + cof_u[i - 1 + N] * basis_shinshkin_1(N, i, tao, 1, x));

	int i = N;
	if ((x > (i - 1.0) * he_1) && (x <= tao))
		return temp + cof_u[i - 1] * (x - (i - 1.0) * he_1) / he_1;
	else if ((x > tao) && (x < tao + he_2))
		return temp + cof_u[i - 1] * (tao + he_2 - x) / he_2;
	else
		return temp;
}

//离散计算无穷范数误差；
double error_max_shinshkin_1(int N, double tao, vector<double> cof_u)
{
	int n = 1000;
	double delta_x = 1.0 / n;
	vector<double> error(n + 1, 0);
	vector<double> u(n - 1, 0);
	for (int i = 0; i <= n; i++)
		error[i] = fabs(u_2(i * delta_x) - uh_basis_shinshkin_1(N, tao, cof_u, i * delta_x));
	double max = error[0];
	for (int i = 1; i <= n; i++)
		if (error[i] > max)
			max = error[i];
	return max;
}
*/
//基函数计算求解，输出无穷范数误差；
void error_out(vector<vector<double>>& error)
{
	for (int i = 0; i < 6; i++)
	{
		//区间离散程度；
		int N = pow(2, i) * 10;
		double he = 1.0 / N;
		double epsilon = 1e-7;
		double tao = 1.0 - 2 * epsilon * log(N);
		//系数矩阵K；
		vector<vector<double>> cof_k_1(N - 1, vector<double>(N - 1, 0));
		make_cof_k(N, cof_k_1);
		//右侧函数f与基内积向量；
		vector<double> cof_f_1(N - 1, 0);
		make_cof_f(N, cof_f_1);
		//离散解的系数向量；
		vector<double> cof_u_1(N - 1, 0);
		gauss(cof_k_1, cof_f_1, cof_u_1);

		//共轭梯度法求解；
		//int times = 0;             //共轭梯度迭代法迭代次数；
		//double iteration_error = 0;//共轭梯度迭代法的输出参考误差；
		//Conjugate_Gradient(N - 1, cof_k_1, cof_f_1, 1e-20, cof_u_1, times, iteration_error);
		//times = 0;           //共轭梯度迭代法迭代次数；
		//iteration_error = 0;//共轭梯度迭代法的输出参考误差；
		//Conjugate_Gradient(2*N - 1, cof_k_shinshkin, cof_f_shinshkin, 1e-10, cof_u_shinshkin, times, iteration_error);

		//计算无穷范数作为误差；
		error[i][0] = error_max_1(N, cof_u_1);
		error[i][2] = error_L1_1(N, cof_u_1);
		//error[i][4] = error_max_shinshkin_1(N, tao, cof_u_shinshkin);
		//error[i][6] = error_L1_shinshkin_1(N, cof_u_2);
		if (i > 0)
		{
			error[i][1] = log(error[i - 1][0] / error[i][0]) / log(2);
			error[i][3] = log(error[i - 1][2] / error[i][2]) / log(2);
			//error[i][5] = log(error[i - 1][4] / error[i][4]) / log(2);
			//error[i][7] = log(error[i - 1][6] / error[i][6]) / log(2);
		}
	}
}