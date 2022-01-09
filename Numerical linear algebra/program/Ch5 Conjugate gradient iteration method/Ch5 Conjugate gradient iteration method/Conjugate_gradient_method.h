#include <iostream>
#include<iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <ctime>

using namespace std;

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
			cout << matrix[i][j] << "   ";
		cout << endl;
	}
}

//生成Hilbert方阵；
void Hilbert_matrix(int n, vector<vector<double>>& matrix)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			matrix[i][j] = 1.0 / (i + j + 1.0);
}

//Ch5.1所需系数矩阵以及向量b的生成；
void Ch5_1(int dim, vector<vector<double>>& A, vector<double>& b)
{
	
	int i, j, k;
	double h = 1.0 / 20;
	for (i = 0; i < 361; i++) A[i][i] = 1.0 + h * h / 4.0;
	for (k = 0; k < 19; k++)
	{
		for (i = 19 * k + 1; i < 19 * k + 18; i++)
		{
			A[i][i + 1] = -0.25;
			A[i][i - 1] = -0.25;
		}
		A[19 * k][19 * k + 1] = -0.25, A[19 * k + 18][19 * k + 17] = -0.25;
	}
	for (k = 1; k < 18; k++) for (i = 19 * k; i < 19 * k + 19; i++) A[i][i + 19] = -0.25, A[i][i - 19] = -0.25;
	for (i = 0; i < 19; i++) A[i][i + 19] = -0.25;
	for (i = 342; i < 361; i++) A[i][i - 19] = -0.25;

	for (i = 1; i < 20; i++) for (j = 1; j < 20; j++) b[19 * (i - 1) + j - 1] = h * h * (sin(h * h * i * j)) / 4.0;
	for (j = 1; j < 20; j++) b[j - 1] += j * j / 1600.0, b[341 + j] += 0.25 + j * j / 1600.0;//i=1,19    19(i-1)+j-1
	for (i = 1; i < 20; i++) b[19 * i - 19] += i * i / 1600.0, b[19 * i - 1] += 0.25 + i * i / 1600.0;//j=1,19

	
	/*int n = (dim - 1) * (dim - 1);
	double h = 1.0 / dim;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			matrix[i][i] = 1.0 + h * h / 4.0;
			if (j == i + 1 || j == i + dim - 1 || j + 1 == i || j + dim - 1 == i) matrix[i][j] = -0.25;
			else matrix[i][j] = 0;
		}
	}
	for (int i = 1; i < dim - 1; i++)//处理不该赋0的位置；
	{
		matrix[i * (dim - 1) - 1][i * (dim - 1)] = 0;
		matrix[i * (dim - 1)][i * (dim - 1) - 1] = 0;
	}
	for (int j = 1; j < dim; j++)//初始b；
		for (int i = 1; i < dim; i++)
			b[(dim-1)*(j-1)+i-1] = h * h * sin(i * j * h * h) / 4.0;
	for (int i = 1; i < dim; i++)//调整b；
	{
		b[i-1] += (i * i * h * h);//书中网格下侧边界值处理；
		b[(i - 1) * (dim - 1)] += (i * i * h * h);//左侧边界值处理；
		b[i * (dim - 1) - 1] += (dim - 1) * (dim - 1) * h * h + i * i * h * h;//右侧边界值处理；
		b[(dim-1)*(dim-2)+i-1]+= (dim - 1) * (dim - 1) * h * h + i * i * h * h;//上侧边界值处理；
	}
	*/
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
void Vector_subtraction(vector<double> x, vector<double> y, vector<double>& dx)
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
	vector<double> x0(dim, 1);
	vector<double> r0(dim, 0), r(dim, 0), p0(dim, 0);
	Matrix_times_vector(matrix, x0, r0);
	Vector_subtraction(b, r0, r0);//顺序？
	r = r0;
	x = x0;
	times = 0;
	while (vector_2_norm(r) > e && times<10000)
	{
		times++;
		if (times == 1) p0 = r0;
		else 
		{
			double beta = Inner_product_of_vectors(r, r) / Inner_product_of_vectors(r0, r0);
			xk(p0, dim, r, beta, p0);
		}
		double alpha= Inner_product_of_vectors(r, r) / pAp(matrix,p0);
		xk(x, dim, x, alpha, p0);
		r0 = r;//储存rk-1;
		rk(r, dim, r0, p0, alpha, matrix);
	}
	error = vector_2_norm(r);
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

//Jacobi迭代矩阵生成，并处理b；
void Jacobi_matrix(int dim, vector<vector<double>>& matrix, vector<double>& value)
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
void Jacobi_iteration(int dim, vector<vector<double>>& matrix, vector<double>& x0, vector<double> value, double error, int& times)
{
	vector<double> x(dim), dx(dim, 0.0), b(dim);
	b = value;//记录b
	x = x0;  //迭代所需的过渡向量；
	times = 0;
	Jacobi_matrix(dim, matrix, b);//生成Jacobi迭代矩阵，并处理b；
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
