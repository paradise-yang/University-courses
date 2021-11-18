#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <time.h>

using namespace std;

//生成n维Hilbert方阵；
void Hilbert_matrix(int n, vector<vector<double>>& matrix)
{
	int i = 0, j = 0;
	vector<double> row(n);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			row[j] = 1.0 / (i + j + 1.0);
		matrix.push_back(row);
	}
}

//生成n维题干系数矩阵；
void creat_matrix(int dim, vector<vector<double>>& matrix)
{
	int i = 0, j = 0;
	vector<double> row(dim);
	for (i = 0; i < dim; i++)
	{
		for (j = 0; j < dim - 1; j++)
		{
			if (i > j) row[j] = -1.0;
			else if (i == j) row[j] = 1.0;
			else row[j] = 0;
		}
		row[dim - 1] = 1.0;
		matrix.push_back(row);
	}
}

//随机生成n维向量；
void Vector_creat_random(int dim, vector<double>& row)
{
	for (int i = 0; i < dim; i++)
	{
		//srand(i);
		//row[i] = sqrt(2)*rand();
		row[i] = sqrt(2) *(i+1.0);
	}
}

//向量输出；
void Vector_cout(vector<double> x)
{
	int dim;
	dim = size(x);
	for (int i = 0; i < dim; i++)
		cout << x[i] << " ";
	cout << endl;
}

//矩阵输出；
void Matrix_cout(vector<vector<double>> matrix)
{
	int i = 0, j = 0, dim;
	dim = size(matrix);
	for (i = 0; i < dim; i++)
	{
		for (j = 0; j < dim; j++)
			cout << matrix[i][j] << " ";
		cout << endl;
	}
}

//矩阵转置；
void Matrix_transpose(vector<vector<double>> matrix, vector<vector<double>>& matrixT)
{
	int dim, i, j;
	dim = size(matrix);
	vector<double> row(dim);
	for (j = 0; j < dim; j++)
	{
		for (i = 0; i < dim; i++)
			row[i] = matrix[i][j];
		matrixT.push_back(row);
	}
}

//矩阵乘以列向量；
void Matrix_times_vector(vector<vector<double>> matrix, vector<double> column, vector<double>& value)
{
	int dim, i, j;
	dim = size(matrix);
	for (i = 0; i < dim; i++)
	{
		value[i] = 0;
		for (j = 0; j < dim; j++)
			value[i] += matrix[i][j] * column[j];
	}
}

//向量之间作减法；
void Vector_subtraction(vector<double> &x, vector<double> y)
{
	int dim = 0;
	dim = size(x);
	for (int i = 0; i < dim; i++)
		x[i] = x[i] - y[i];
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

//生成向量的符号向量；
void Vector_sign(vector<double> column, vector<double>& sign)
{
	int n, i;
	n = size(column);
	for (i = 0; i < n; i++)
	{
		if (column[i] > 0) sign[i] = 1.0;
		else if (column[i] < 0) sign[i] = -1.0;
		else sign[i] = 0;
	}
}

//生成第k个基向量ek；
void base_vector(int n,vector<double>& Ek, int k)
{
	for (int i = 0; i < n; i++)
	{
		if (i == k) Ek[i] = 1.0;
		else Ek[i] = 0;
	}
}

//向量的无穷范数，并记录所取到位置（实际位置-1）；
void Vector_infinity_norm(vector<double> x, double& norm, int& j)
{
	int i, n;
	norm = fabs(x[0]);
	n = size(x);
	j = 0;
	for (i = 1; i < n; i++)
	{
		if (fabs(x[i]) >= norm)
		{
			norm = fabs(x[i]);
			j = i;
		}
	}
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

//向量的1-范数；
double Vector_1_norm(vector<double> x)
{
	int i, n;
	double temp = 0;
	n = size(x);
	for (i = 0; i < n; i++)
		temp += fabs(x[i]);
	return temp;
}

//计算矩阵的1-范数，并返回范数值；
double Matrix_1_norm(vector<vector<double>> matrix)
{
	int maxj = 0, i, dim;
	double temp = 0, norm = 0;
	dim = size(matrix);                       //记录矩阵（方阵）维数；
	vector<double> w(dim), x(dim), v(dim), z(dim);
	vector<vector<double>> matrixT;
	Matrix_transpose(matrix, matrixT);        //记录矩阵的转置；
	for (i = 0; i < dim; i++)                 //初始化x；
		x[i] = 1.0 / dim;
	while (1)
	{
		Matrix_times_vector(matrix, x, w);    //初始化w；
		Vector_sign(w, v);                    //记录符号向量；
		Matrix_times_vector(matrixT, v, z);   //计算z；
		temp = Inner_product_of_vectors(z, x);//计算z与x内积；
		Vector_infinity_norm(z, norm, maxj);  //计算z的无穷范数，并记录所取位置；
		if (norm <= temp)
			return Vector_1_norm(w);          //返回局部最大值w的1-范数；
		else base_vector(dim, x, maxj);       //取无穷范数所取到的列数生成基向量；
	}
}

//列主元Guass消去求解LU，并同时置换；
int Guass_Line(int n, vector<vector<double>>& matrix, vector<double>& value)
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

//前代法解单位下三角方程组；
int lower_unit_triangular(int n, vector<vector<double>> matrix, vector<double>& value)
{
	int i, j;
	for (i = 0; i < n - 1; i++)
		for (j = i + 1; j < n; j++)
			value[j] = value[j] - value[i] * matrix[j][i];
	return 0;
}

//回代法解上三角方程组；
int upper_triangular(int n, vector<vector<double>> matrix, vector<double>& value)
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

//列主元的Guass消去法解方程，并且不改变系数矩阵，将解储存在value中；
void Guass_Line_solution(vector<vector<double>> matrix, vector<double> column, vector<double>& value)
{
	int i, dim;
	dim = size(matrix);
	vector<vector<double>> matrix_temp;
	for (i = 0; i < dim; i++)
	{
		matrix_temp.push_back(matrix[i]);
		value[i] = column[i];
	}                                               //生成暂替系数矩阵；
	Guass_Line(dim, matrix_temp, value);           //列主元Guass消去；
	lower_unit_triangular(dim, matrix_temp, value);//解单位下三角得出y;
	upper_triangular(dim, matrix_temp, value);     //解上三角得出x;
}

//计算矩阵的逆的1-范数，并返回范数值；
double Matrix_inverse_1_norm(vector<vector<double>> matrix)
{
	int maxj = 0, i, dim;
	double temp = 0, norm = 0;
	dim = size(matrix);                           //记录矩阵（方阵）维数；
	vector<double> w(dim), x(dim), v(dim), z(dim);
	vector<vector<double>> matrixT;
	Matrix_transpose(matrix, matrixT);            //记录矩阵的转置；
	double n = 0;
	for (i = 0; i < dim; i++)                     //初始化x；
		x[i] = 1.0 / dim;
	while (1)
	{
		n = n + 1;
		Guass_Line_solution(matrix, x, w);        //初始化w；
		Vector_sign(w, v);                        //记录符号向量；
		Guass_Line_solution(matrixT, v, z);       //计算z；
		temp = Inner_product_of_vectors(z, x);    //计算z与x内积；
		Vector_infinity_norm(z, norm, maxj);      //计算z的无穷范数，并记录所取位置；
		if (norm <= temp) return Vector_1_norm(w);//返回局部最大值w的1-范数；
		else base_vector(dim, x, maxj);           //取无穷范数所取到的列数生成基向量；
	}
}

//计算并返回矩阵的无穷范数的条件数；
double Conditional_number_infinity(vector<vector<double>> matrix)
{
	int dim;
	dim = size(matrix);
	vector<vector<double>> matrixT;
	Matrix_transpose(matrix, matrixT);
	return Matrix_1_norm(matrixT) * Matrix_inverse_1_norm(matrixT);
}