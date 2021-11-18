#include <iostream>
#include<iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <ctime>

using namespace std;

//生成带状系数矩阵；
void band_matrix(int n, vector<vector<double>>& matrix, double epsilon, double h) 
{
	double eh = epsilon + h, eeh = epsilon + epsilon + h;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j) matrix[i][j] = 0 - eeh;
			else if (j == i + 1) matrix[i][j] = eh;
			else if (i == j + 1) matrix[i][j] = epsilon;
			else matrix[i][j] = 0;
		}
	}
}

//Ch4.1 精确函数，得到精确函数解yi；
void F(int dim, double epsilon, double a, vector<double>& y)
{
	double h = -1 / (dim * epsilon);
	double c = (1 - a) / (1 - exp(-1 / epsilon));
	for (int i = 0; i < dim - 1; i++)
		y[i] = c * (1 - exp((i + 1) * h)) + a * (i + 1) / dim;
}

//Ch4.2所需函数f；
double f(int N, int i, int j)
{
	return (i + j) / N;
}

//Ch4.2所需函数g；
double g(int N, int i, int j)
{
	return exp(i * j / (N * N));
}

//Ch4.2相关矩阵生成，U为初始迭代矩阵，为n+1阶，边界储存为1；
void Ch4_2(int dim, vector<vector<double>>& U, vector<vector<double>>& F, vector<vector<double>>& G)
{
	double h = dim ^ 2;
	for (int i = 0; i < dim + 1; i++)
	{
		for (int j = 0; j < dim + 1; j++)
		{
			F[i][j] = f(dim, i, j) / h;//储存的是h^2*f(ij)；
			G[i][j] = g(dim, i, j) / h;//储存的是h^2*g(ij)；
			if (j == 0 || i == 0 || i == dim || j == dim) U[i][j] = 1.0;
			else U[i][j] = 1.0 / dim;
		}
	}
}

//向量输出；
void Vector_cout(vector<double> x)
{
	int dim = size(x);
	for (int i = 0; i < dim; i++)
		cout << setprecision(4) << x[i] << "   ";
	cout << endl;
}

//矩阵输出；
void Matrix_cout(vector<vector<double>> matrix)
{
	int dim = size(matrix);
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
			cout << matrix[i][j] << " ";
		cout << endl;
	}
}

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
void Vector_subtraction(vector<double> x, vector<double> y, vector<double> &dx)
{
	int dim = 0;
	dim = size(x);
	for (int i = 0; i < dim; i++)
		dx[i] = x[i] - y[i];
}

//向量的无穷范数，并返回范数值；
double vector_infinity_norm(vector<double> x)
{
	int i, n;
	double norm = fabs(x[0]);
	n = size(x);
	for (i = 1; i < n; i++)
	{
		if (fabs(x[i]) > norm)
			norm = fabs(x[i]);
	}
	return norm;
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

//矩阵乘法,并储存在第二个矩阵里面；
void Matrix_times(vector<vector<double>> house, vector<vector<double>>& matrix, int m, int n)
{
	vector<vector<double>> temp(m, vector<double>(n, 0));
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < m; k++)
				temp[i][j] += house[i][k] * matrix[k][j];
		}
	}
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			matrix[i][j] = temp[i][j];
}

//Jacobi迭代矩阵生成，并处理b；
void Jacobi_matrix(int dim, vector<vector<double>>& matrix, vector<double> &value)
{
	for (int i = 0; i < dim; i++)
	{
		double temp = matrix[i][i];
		matrix[i][i] = 0;
		value[i] = value[i] / temp;
		for (int j = 0; j < dim; j++)
			matrix[i][j] = -1 * matrix[i][j] / temp;
	}
}

//Jacobi迭代法,x0为初始向量，value为右端b，error为最大误差（有效数字）；
void Jacobi_iteration(int dim, vector<vector<double>>& matrix, vector<double> &x0, vector<double> value, double error, int &times)
{
	vector<double> x(dim), dx(dim, 0.0), b(dim);
	b = value;//记录b
	x = x0;  //迭代所需的过渡向量；
	times = 0;
	Jacobi_matrix(dim, matrix,b);//生成Jacobi迭代矩阵，并处理b；
	do {
		times++;
		Matrix_times_vector(matrix, x0, x);
		for (int i = 0; i < dim; i++)
			x[i] += b[i];
		Vector_subtraction(x0, x, dx);
		x0 = x;
	} while (vector_infinity_norm(dx) >= error);
}

//Guass-Seidel迭代法,x0为初始向量，value为右端b，error为最大误差（有效数字）；
void Gauss_Seidel_iteration(int dim, vector<vector<double>>& matrix, vector<double>& x0, vector<double> value, double error, int& times)
{
	vector<double> x(dim), dx(dim, 0.0), b(dim);
	x = x0;
	b = value;
	Jacobi_matrix(dim, matrix, b);
	times = 0;
	do {
		times++;
		for (int i = 0; i < dim; i++)
		{
			double temp = 0.0;
			for (int j = 0; j < dim; j++)
				temp += matrix[i][j] * x0[j];
			x0[i] = temp + b[i];
		}
		Vector_subtraction(x0, x, dx);
		x = x0;
	} while (vector_infinity_norm(dx) >= error);
}

//处理SOR所需的迭代矩阵以及b；
void SOR_matrix(int dim, double w, vector<vector<double>>& matrix, vector<double>& value)
{
	for (int i = 0; i < dim; i++)
	{
		double temp = w / matrix[i][i];
		value[i] *= temp;
		for (int j = 0; j < dim; j++)
			matrix[i][j] = -1 * matrix[i][j] * temp;
		matrix[i][i] = 1.0 - w;
	}
}

//SOR超松弛迭代法,x0为初始向量，value为右端b，error为最大误差（有效数字）；
void SOR_iteration(int dim, double w, vector<vector<double>>& matrix, vector<double>& x0, vector<double> value, double error, int& times)
{
	vector<double> x(dim), dx(dim, 0.0), b(dim);
	x = x0;
	b = value;
	SOR_matrix(dim, w, matrix, b);
	times = 0;
	do {
		times++;
		for (int i = 0; i < dim; i++)
		{
			double temp = 0.0;
			for (int j = 0; j < dim; j++)
				temp += matrix[i][j] * x0[j];
			x0[i] = temp + b[i];
		}
		Vector_subtraction(x0, x, dx);
		x = x0;
	} while (vector_infinity_norm(dx) >= error);
}

//矩阵作减法；
void matrix_subtraction(vector<vector<double>> M, vector<vector<double>> m, vector<vector<double>>& dm)
{
	int dim = size(M);
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			dm[i][j] = M[i][j] - m[i][j];
}

//矩阵作为向量求2范数，并返回范数值；
double matrix_as_vector_2_norm(vector<vector<double>> matrix)
{
	int dim = size(matrix);
	double norm = 0;
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			norm = norm + matrix[i][j] * matrix[i][j];
	return sqrt(norm);
}

//Jacobi迭代求解偏微分方程离散解；
void Jacobi_PDE(int dim, int& times, double error)
{
	vector<vector<double>> U(dim + 1, vector<double>(dim + 1, 0)), U0(dim + 1, vector<double>(dim + 1, 0)), dU(dim + 1, vector<double>(dim + 1, 0));
	vector<vector<double>> F(dim + 1, vector<double>(dim + 1, 0)), G(dim + 1, vector<double>(dim + 1, 0));
	Ch4_2(dim, U, F, G);
	U0 = U;
	times = 0;
	do {
		times++;
		for (int j = 1; j < dim; j++)
			for (int i = 1; i < dim; i++)
				U0[i][j] = (U[i][j - 1] + U[i][j + 1] + U[i - 1][j] + U[i + 1][j] + F[i][j]) / (4 + G[i][j]);
		matrix_subtraction(U, U0, dU);
		U = U0;
	} while (matrix_as_vector_2_norm(dU) >= error);
}

//Guass-Seidel迭代求解偏微分方程离散解；
void Gauss_Seidel_PDE(int dim, int &times, double error)
{
	vector<vector<double>> U(dim + 1, vector<double>(dim + 1, 0)), U0(dim + 1, vector<double>(dim + 1, 0)), dU(dim + 1, vector<double>(dim + 1, 0));
	vector<vector<double>> F(dim + 1, vector<double>(dim + 1, 0)), G(dim + 1, vector<double>(dim + 1, 0));
	Ch4_2(dim, U, F, G);
	U0 = U;
	times = 0;
	do {
		times++;
		for (int j = 1; j < dim; j++)
			for (int i = 1; i < dim; i++)
				U[i][j] = (U[i][j - 1] + U[i][j + 1] + U[i - 1][j] + U[i + 1][j] + F[i][j]) / (4 + G[i][j]);
		matrix_subtraction(U, U0, dU);
		U0 = U;
	} while (matrix_as_vector_2_norm(dU)>=error);
}

//SOR迭代求解偏微分方程离散解；
void SOR_PDE(int dim, int& times, double error, double w)
{
	vector<vector<double>> U(dim + 1, vector<double>(dim + 1, 0)), U0(dim + 1, vector<double>(dim + 1, 0)), dU(dim + 1, vector<double>(dim + 1, 0));
	vector<vector<double>> F(dim + 1, vector<double>(dim + 1, 0)), G(dim + 1, vector<double>(dim + 1, 0));
	Ch4_2(dim, U, F, G);
	U0 = U;
	times = 0;
	do {
		times++;
		for (int j = 1; j < dim; j++)
		{
			for (int i = 1; i < dim; i++)
				U[i][j] = (1 - w) * U[i][j] + w * (U[i][j - 1] + U[i][j + 1] + U[i - 1][j] + U[i + 1][j] + F[i][j]) / (4 + G[i][j]);
		}
		matrix_subtraction(U, U0, dU);
		U0 = U;
	} while (matrix_as_vector_2_norm(dU) >= error);
}