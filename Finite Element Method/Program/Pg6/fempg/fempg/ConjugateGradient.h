#pragma once
#include<iostream>
#include <iomanip>
#include<vector>
#include<cmath>
#include <chrono>

using namespace std;
#define ITERATION_MAX 1000
#define PRECISION 1e-8

//顺哥的///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef unsigned index;
typedef std::vector<std::vector<double>> matrix;

//相关运算
//向量点乘向量
double vectorMultiplyvector(const vector<double> &u, const vector<double> &v)
{
	index n1 = u.size();
	index n2 = v.size();
	double s = 0;
	if (n1 != n2) {
		printf("vectorMultiplyvector:向量点乘尺寸异常\n");
		return s;
	}
	for (index i = 0; i < n1; i++) {
		s = s + u[i] * v[i];
	}
	return s;
}
//向量加法重载
vector<double> operator+(const vector<double> &u, const vector<double> &v)
{
	index n1 = u.size();
	index n2 = v.size();
	if (n1 != n2) {
		printf("operator+:向量相加尺寸不同\n");
		return u;
	}

	vector<double> w(n1);
	for (index i = 0; i < n1; i++)
		w[i] = u[i] + v[i];

	return w;
}
//向量减法重载
vector<double> operator-(const vector<double> &u, const vector<double> &v)
{
	index n1 = u.size();
	index n2 = v.size();
	if (n1 != n2) {
		printf("operator-:向量相减尺寸不同\n");
		return u;
	}

	vector<double> w(n1);
	for (index i = 0; i < n1; i++)
		w[i] = u[i] - v[i];

	return w;
}
//向量数乘重载
vector<double> operator*(double x, const vector<double> &u)
{
	index n = u.size();
	vector<double> w(n);
	for (index i = 0; i < n; i++) {
		w[i] = x * u[i];
	}
	return w;
}
vector<double> quick_multiply_3_line(const matrix &A, const vector<double> &b)
{
	index n = A.size();
	if (n != A[0].size() || n != b.size()) {
		printf("error\n");
		return b;
	}
	vector<double> result(n, 0);

	for (index i = 0; i < n - 1; i++) {
		result[i] += (A[i][i] * b[i]);
		result[i] += (A[i][i + 1] * b[i + 1]);
		result[i + 1] += (A[i + 1][i] * b[i]);
	}
	result[n - 1] += (A[n - 1][n - 1] * b[n - 1]);

	return result;
}
vector<double> Conjugate_Gradient(const matrix &A, const vector<double> &b, const vector<double> &start, double precision, index num_max)
{
	if (A.size() != A[0].size() || A[0].size() != b.size() || b.size() != start.size()) {
		printf("共轭梯度法:矩阵与向量尺寸有问题\n");
		return b;
	}

	vector<double> x = start;
	vector<double> r(b), p(b), w(b), temp(b);

	index k = 0;
	index n = b.size();
	index k_max = 10 * n;

	if (num_max > 0 && k_max < num_max) k_max = num_max; //妥善处理最大迭代次数的参数

	double rou = 0, pre_rou = 0;
	double beta = 0, alpha = 0;
	double temp_T = 0;
	double epsilon = 1e-8;
	if (precision > 0) epsilon = precision; //妥善处理精度参数

	double up = epsilon * epsilon * (vectorMultiplyvector(b, b)); //误差上限平方

	k = 0;
	// temp = A * x;
	temp = quick_multiply_3_line(A, x);
	r = b - temp;
	rou = vectorMultiplyvector(r, r);
	while (rou > up && k < k_max) {
		k = k + 1;
		if (k == 1) {
			p = r;
		}
		else {
			beta = rou / pre_rou;
			p = r + beta * p;
		}
		// w = A * p; //每次迭代只有这里有一次矩阵乘以向量，是最大的计算量
		w = quick_multiply_3_line(A, p);
		temp_T = vectorMultiplyvector(p, w);
		alpha = rou / temp_T;
		x = x + alpha * p;
		r = r - alpha * w;

		pre_rou = rou;
		rou = vectorMultiplyvector(r, r);
	}
	if (k >= k_max) {
		printf("共轭梯度法:超过最大迭代次数%d\n", k_max);
	}

	return x;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//output a matrix
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
//output a vector
void vector_output(vector<double> vec)
{
	//cout << "向量输出如下：" << endl;
	for (int i = 0; i < vec.size(); i++)
		cout << vec[i] << "\t";
	cout << endl;
}
//matrix*vector
vector<double> matrix_times_vector(vector<vector<double>> matrix, vector<double> column)
{
	if (matrix[0].size() != column.size())
		cout << "矩阵的列数与向量维数不相等，无法做乘法！" << endl;
	///*一般矩阵情况；
	int row_number = matrix.size(), column_number = matrix[0].size();
	vector<double> value(row_number, 0);
	for (int i = 0; i < row_number; i++)
		for (int j = 0; j < column_number; j++)
			value[i] += matrix[i][j] * column[j];
	return value;
	//*/
	/*三对角矩阵；
	int row_number = matrix.size();
	vector<double> value(row_number, 0);
	value[0] = matrix[0][0] * column[0] + matrix[0][1] * column[1];
	for (int i = 1; i < row_number - 1; i++)
		for (int j = i - 1; j <= i + 1; j++)
			value[i] += matrix[i][j] * column[j];
	value[row_number - 1] = matrix[row_number - 1][row_number - 2] * column[row_number - 2] + matrix[row_number - 1][row_number - 1] * column[row_number - 1];
	return value;
	//*/
}
//三对角矩阵快速乘法；
vector<double> matrix_times_vector_quick(vector<vector<double>> matrix, vector<double> column)
{
	if (matrix[0].size() != column.size())
		cout << "矩阵的列数与向量维数不相等，无法做乘法！" << endl;
	int row_number = matrix.size();
	vector<double> value(row_number, 0);
	value[0] = matrix[0][0] * column[0] + matrix[0][1] * column[1];
	for (int i = 1; i < row_number - 1; i++)
		for (int j = i - 1; j <= i + 1; j++)
			value[i] += matrix[i][j] * column[j];
	value[row_number - 1] = matrix[row_number - 1][row_number - 2] * column[row_number - 2] + matrix[row_number - 1][row_number - 1] * column[row_number - 1];
	return value;
}

//vector+vector
vector<double> vector_add_vector(vector<double> x, vector<double> y)
{
	if (x.size() != y.size())
		cout << "两个向量维数不相等，无法做加法！" << endl;
	int dim = x.size();
	vector<double> add(dim, 0);
	for (int i = 0; i < dim; i++)
		add[i] = x[i] + y[i];
	return add;
}
//vector-vector
vector<double> vector_subtract_vector(vector<double> x, vector<double> y)
{
	if (x.size() != y.size())
		cout << "两个向量维数不相等，无法做减法！" << endl;
	int dim = x.size();
	vector<double> dx(dim, 0);
	for (int i = 0; i < dim; i++)
		dx[i] = x[i] - y[i];
	return dx;
}
//vector^T*vector
double vector_inner_product(vector<double> x, vector<double> y)
{
	if (x.size() != y.size())
		cout << "两个向量维数不相等，无法做内积！" << endl;
	int dim = x.size();
	double temp = 0;
	for (int i = 0; i < dim; i++)
		temp += x[i] * y[i];
	return temp;
}

/////////以下为CG特定的辅助函数；
//p^t*A*p
//double pAp(vector<vector<double>> matrix, vector<double> column)
//{
//	int dim = matrix.size();
//	return vector_inner_product(column, matrix_times_vector(matrix, column));
//}
//the direction of iteration 
vector<double> xk(vector<double> x0, double alpha, vector<double> p0)
{
	int dim = x0.size();
	vector<double> x(dim, 0);
	for (int i = 0; i < dim; i++)
		x[i] = x0[i] + alpha * p0[i];
	return x;
}
//the residual vector at step k
vector<double> rk(vector<double> r0, vector<double> p0, double a0, vector<vector<double>> matrix)
{
	int dim = r0.size();
	vector<double> temp(dim, 0);
	temp = matrix_times_vector(matrix, p0);
	for (int i = 0; i < dim; i++)
		temp[i] = r0[i] - a0 * temp[i];
	return temp;
}
//Conjugate Gradient method, return the time taken to solve the equation and the solution
vector<double> conjugate_gradient(vector<vector<double>> matrix, vector<double> b, vector<double> start, double precision, int& iteration_times, double& time)
{
	if (matrix.size() != matrix[0].size() || matrix[0].size() != b.size() || start.size() != b.size())
		cout << "矩阵行数与列数不相等，或，矩阵的列数与右侧向量维数不相等，或初始迭代向量维数与右侧向量维数不等，无法求解！" << endl;
	int dim = matrix.size(), iteration = 0;
	if (iteration_times < 10 * dim) iteration_times = 10 * dim;//deal with the number of iteration
	if (iteration_times < ITERATION_MAX) iteration_times = ITERATION_MAX;
	double up_precision = precision * precision*vector_inner_product(b, b);
	vector<double> x = start, p0(b), w(b);
	vector<double> r0 = matrix_times_vector(matrix, x);
	vector<double> r = b - r0;//vector_subtract_vector(b, r0);
	double alpha = 0, beta = 0, rou = vector_inner_product(r, r);
	auto start_time = std::chrono::steady_clock::now();
	while (rou > up_precision&&iteration < iteration_times)
	{
		iteration++;
		if (iteration == 1) p0 = r;
		else beta = rou / vector_inner_product(r0, r0), p0 = r + beta * p0;//xk(r, beta, p0);
		w = matrix_times_vector(matrix, p0);
		alpha = rou / vector_inner_product(p0, w);//alpha = rou / pAp(matrix, p0);
		x = x + alpha * p0;//xk(x, alpha, p0);
		r0 = r;
		r = r0 - alpha * w;//rk(r0, p0, alpha, matrix);r - alpha * w
		rou = vector_inner_product(r, r);
	}
	auto end_time = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end_time - start_time;
	time = elapsed_seconds.count();
	iteration_times = iteration;
	return x;
}
//Conjugate Gradient method,三对角；
vector<double> conjugate_gradient_quick(vector<vector<double>> matrix, vector<double> b, vector<double> start, double precision, int& iteration_times, double& time)
{
	if (matrix.size() != matrix[0].size() || matrix[0].size() != b.size() || start.size() != b.size())
		cout << "矩阵行数与列数不相等，或，矩阵的列数与右侧向量维数不相等，或初始迭代向量维数与右侧向量维数不等，无法求解！" << endl;
	int dim = matrix.size(), iteration = 0;
	if (iteration_times < 10 * dim) iteration_times = 10 * dim;//deal with the number of iteration
	if (iteration_times < ITERATION_MAX) iteration_times = ITERATION_MAX;
	double up_precision = precision * precision*vector_inner_product(b, b);
	vector<double> x = start, p0(b), w(b);
	vector<double> r0 = matrix_times_vector_quick(matrix, x);
	vector<double> r = b - r0;//vector_subtract_vector(b, r0);
	double alpha = 0, beta = 0, rou = vector_inner_product(r, r);
	auto start_time = std::chrono::steady_clock::now();
	while (rou > up_precision&&iteration < iteration_times)
	{
		iteration++;
		if (iteration == 1) p0 = r;
		else beta = rou / vector_inner_product(r0, r0), p0 = r + beta * p0;//xk(r, beta, p0);
		w = matrix_times_vector_quick(matrix, p0);
		alpha = rou / vector_inner_product(p0, w);//alpha = rou / pAp(matrix, p0);
		x = x + alpha * p0;//xk(x, alpha, p0);
		r0 = r;
		r = r0 - alpha * w;//rk(r0, p0, alpha, matrix);r - alpha * w
		rou = vector_inner_product(r, r);
	}
	auto end_time = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end_time - start_time;
	time = elapsed_seconds.count();
	iteration_times = iteration;
	return x;
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

int lower_unit_triangular(int n, vector<vector<double>> matrix, vector<double>& value)  //单位下三角方阵前代法求解
{
	int i, j;
	for (i = 0; i < n - 1; i++)
		for (j = i + 1; j < n; j++)
			value[j] = value[j] - value[i] * matrix[j][i];
	return 0;
}
int upper_unit_triangular(int n, vector<vector<double>> matrix, vector<double>& value)  //单位上三角方阵回代法求解
{
	int i, j;
	for (i = n - 1; i > 0; i--)
		for (j = 0; j < i; j++)
			value[j] = value[j] - value[i] * matrix[j][i];
	return 0;
}
int D(int n, vector<vector<double>> matrix, vector<double>& value)
{
	int i = 0;
	for (i = 0; i < n; i++)
		value[i] = value[i] / matrix[i][i];
	return 0;
}
int Cholesky_improve(int n, vector<vector<double>>& matrix)//改进平方根法计算Cholesky分解；
{
	int i = 0, j = 0, k = 0;
	vector<double> v(n);
	double temp = 0;
	for (j = 0; j < n; j++)
	{
		for (i = 0; i < j; i++)
			v[i] = matrix[j][i] * matrix[i][i];
		for (i = 0; i < j; i++)
			matrix[j][j] = matrix[j][j] - matrix[j][i] * v[i];
		for (k = j + 1; k < n; k++)
		{
			temp = 0;
			for (i = 0; i < j; i++)
				temp += matrix[k][i] * v[i];
			matrix[k][j] = (matrix[k][j] - temp) / matrix[j][j];
		}
	}
	for (i = 0; i < n; i++)//将方阵上侧变成L;
		for (j = i + 1; j < n; j++)
			matrix[i][j] = matrix[j][i];
	return 0;
}
vector<double> cholesky(vector<vector<double>> matrix, vector<double> value)
{
	int dim = matrix.size();
	Cholesky_improve(dim, matrix);
	lower_unit_triangular(dim, matrix, value);
	D(dim, matrix, value);
	upper_unit_triangular(dim, matrix, value);
	return value;
}

int lower_triangular(int n, vector<vector<double>> matrix, vector<double>& value)  //（单位）下三角方阵前代法求解
{
	int i, j;
	for (i = 0; i < n - 1; i++)
		for (j = i + 1; j < n; j++)
			value[j] = value[j] - value[i] * matrix[j][i];
	return 0;
}
int upper_triangular(int n, vector<vector<double>> matrix, vector<double>& value)  //上三角方阵回代法求解
{
	int i, j;
	for (i = n - 1; i > 0; i--)
	{
		value[i] = value[i] / matrix[i][i];
		for (j = 0; j < i; j++)
			value[j] = value[j] - value[i] * matrix[j][i];
	}
	value[0] = value[0] / matrix[0][0];
	return 0;
}
int Guass_LU(int n, vector<vector<double>>& matrix)  //未选主元的Guass消去求解LU分解
{
	int k, i, j;
	for (k = 0; k < n - 1; k++)
	{
		for (i = k + 1; i < n; i++)
			matrix[i][k] = matrix[i][k] / matrix[k][k];
		for (i = k + 1; i < n; i++)
			for (j = k + 1; j < n; j++)
				matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
	}
	return 0;
}
int Guass_Line(int n, vector<vector<double>>& matrix, vector<double>& value)  //列主元Guass消去求解LU，并同时置换
{
	int k, i, j, p = 0;
	double max, temp;
	for (k = 0; k < n - 1; k++)
	{
		max = matrix[k][k];//寻找列主元
		p = k;
		for (i = k + 1; i < n; i++)
		{
			if (fabs(matrix[i][k]) > fabs(max))
			{
				max = matrix[i][k];
				p = i;
			}
		}
		if (p != k) //置换；
		{
			vector<double> pc(matrix[k]);
			matrix[k] = matrix[p];
			matrix[p] = pc;
			temp = value[k];
			value[k] = value[p];
			value[p] = temp;
		}

		if (abs(matrix[k][k]) != 0)
		{
			for (i = k + 1; i < n; i++)
				matrix[i][k] = matrix[i][k] / matrix[k][k];
			for (j = k + 1; j < n; j++)
				for (i = k + 1; i < n; i++)
					matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
		}
	}
	return 0;
}
vector<double> gauss_line(vector<vector<double>>& matrix, vector<double>& value)
{
	int dim = matrix.size();
	Guass_Line(dim, matrix, value);      //列主元Guass消去；
	//Guass_LU(dim, matrix);
	lower_triangular(dim, matrix, value);//解下三角得出y;
	upper_triangular(dim, matrix, value);//解上三角得出x;
	return value;
}