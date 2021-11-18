#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

//生成Ch1.1的带状系数方阵；
void band_matrix1(int n, vector<vector<double>>& matrix) 
{
	int i, j;
	vector<double> row(n);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i == j) row[j] = 6;
			else if (j == i + 1) row[j] = 1;
			else if (j == i - 1) row[j] = 8;
			else row[j] = 0;
		}
		matrix.push_back(row);
	}
}

//生成Ch1.2.(1)的带状系数方阵；
void band_matrix2(int n, vector<vector<double>>& matrix)      //生成带状系数方阵
{
	int i, j;
	vector<double> row(n);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i == j) row[j] = 10;
			else if (j == i + 1 || j == i - 1) row[j] = 1;
			else row[j] = 0;
		}
		matrix.push_back(row);
	}
}

//生成n阶Hilbert方阵；
void Hilbert_matrix(int n, vector<vector<double>>& matrix)  //生成Hilbert方阵；
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

//生成Ch3.2所需的系数方阵以及向量b；
void Ch3_2_make(vector<vector<double>>& matrix, vector<double>& b)
{
	vector<double> t(7), row(3);
	t[0] = -1.0; t[1] = -0.75; t[2] = -0.5; t[3] = 0; t[4] = 0.25; t[5] = 0.5; t[6] = 0.75;
	b[0] = 1.0; b[1] = 0.8125; b[2] = 0.75; b[3] = 1.0; b[4] = 1.3125; b[5] = 1.75; b[6] = 2.3125;
	for (int i = 0; i < 7; i++)
	{
		row[0] = t[i]*t[i];
		row[1] = t[i];
		row[2] = 1;
		matrix.push_back(row);
	}
}

//生成Ch3.3所需的系数方阵以及向量b；
void Ch3_3_make(vector<vector<double>>& matrix, vector<double>& b)
{
	vector<vector<double>> A =
	{ {1,4.9176, 1, 3.472, 0.998, 1, 7, 4, 42, 3, 1, 0},
	{1,5.0208, 1, 3.531, 1.5, 2, 7, 4, 62, 1, 1, 0},
	{1,4.5429, 1, 2.275, 1.175, 1, 6, 3, 40,  2, 1, 0},
	{1,4.5573, 1, 4.05, 1.232, 1, 6, 3, 54, 4, 1, 0},
	{1,5.0597, 1, 4.455, 1.121, 1, 6, 3, 42, 3, 1, 0},
	{1,3.891, 1, 4.455, 0.988, 1, 6, 3, 56, 2, 1, 0},
	{1,5.898, 1, 5.85, 1.24, 1, 7, 3, 51, 2, 1,  1},
	{1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
	{1,15.4202, 2.5,  9.8, 3.42, 2, 10, 5, 42, 2, 1, 1},
	{1,14.4598, 2.5, 12.8, 3, 2, 9, 5, 14, 4, 1, 1},
	{1,5.8282, 1, 6.435, 1.225, 2, 6, 3, 32, 1, 1, 0},
	{1,5.3003, 1, 4.9883, 1.552, 1, 6, 3, 30, 1, 2, 0},
	{1,6.2712, 1, 5.52, 0.975, 1, 5, 2, 30, 1, 2, 0},
	{1,5.9592, 1, 6.666, 1.121, 2, 6, 3, 32, 2, 1, 0},
	{1,5.05, 1, 5, 1.02, 0, 5, 2, 46, 4, 1, 1},
	{1,5.6039, 1, 9.52, 1.501, 0, 6, 3, 32, 1, 1, 0},
	{1,8.2464, 1.5, 5.15, 1.664, 2, 8, 4, 50, 4, 1, 0},
	{1,6.6969, 1.5, 6.092, 1.488, 1.5, 7, 3, 22, 1, 1, 1},
	{1,7.7841, 1.5, 7.102, 1.376, 1, 6, 3, 17, 2, 1, 0},
	{1,9.0384, 1, 7.8, 1.5, 1.5, 7, 3, 23, 3, 3, 0},
	{1,5.9894, 1, 5.52, 1.256, 2, 6, 3, 40, 4, 1, 1},
	{1,7.5422, 1.5, 4, 1.69, 1, 6, 3, 22, 1, 1, 0},
	{1,8.7951, 1.5, 9.89, 1.82, 2, 8, 4, 50, 1, 1, 1},
	{1,6.0931, 1.5, 6.7265, 1.652, 1, 6, 3, 44, 4, 1, 0},
	{1,8.3607, 1.5, 9.15, 1.777, 2., 8, 4, 48, 1, 1, 1},
	{1,8.14, 1, 8, 1.504, 2, 7, 3, 3, 1, 3, 0},
	{1,9.1416, 1.5, 7.3262, 1.831, 1.5, 8, 4, 31, 4, 1, 0},
	{1,12, 1.5, 5, 1.2, 2, 6, 3, 30, 3, 1, 1} };
	matrix = A;
	b = { 25.9, 29.5, 27.9, 25.9, 29.9, 29.9, 30.9,
	28.9, 84.9, 82.9, 35.9, 31.5, 31.0, 30.9,
	30.0, 28.9, 36.9, 41.9, 40.5, 43.9, 37.5,
	37.9, 44.5, 37.9, 38.9, 36.9, 45.8, 41.0 };
}

//向量输出；
void Vector_cout(int dim, vector<double> x)
{
	for (int i = 0; i < dim; i++)
		cout << x[i] << " ";
	cout << endl;
}

//矩阵输出；
void Matrix_cout(vector<vector<double>> matrix, int m, int n)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			cout << matrix[i][j] << " ";
		cout << endl;
	}
}

//上三角方阵回代法求解；
void upper_triangular(int n, vector<vector<double>> matrix, vector<double>& value)
{
	int i, j;
	for (i = n - 1; i > 0; i--)
	{
		value[i] = value[i] / matrix[i][i];
		for (j = 0; j < i; j++)
			value[j] = value[j] - value[i] * matrix[j][i];
	}
	value[0] = value[0] / matrix[0][0];
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

//计算Householder变换，使得向量变为e1，返回生成向量v，系数b；
void Householder(vector<double> x, vector<double>& v, double& b)
{
	int dim = 0;
	double max = 0, m = 0, a = 0;
	dim = size(x);
	max = vector_infinity_norm(x);
	for (int i = 1; i < dim; i++)
	{
		v[i] = x[i] / max;
		m += pow(v[i], 2);
	}
	if (m == 0) b = 0;
	else 
	{
		a = sqrt(pow(x[0] / max, 2) + m);
		if (x[0] / max <= 0) v[0] = x[0] / max - a;
		else v[0] = -m / (x[0] / max + a);
		b = 2 * pow(v[0], 2) / (m + pow(v[0], 2));
		for (int i = 1; i < dim; i++)
			v[i] = v[i] / v[0];
		v[0] = 1.0;
	}
}

/*//辅助函数，返回Householder方阵作乘法时某位置所需的元素值；
double return_value(vector<double> v, double b, int i, int k)
{
	if (i == k) return (1.0 - b * v[i] * v[k]);
	else return (0 - b * v[i] * v[k]);
}*/

//Householder方阵生成；
void Householder_make(int n, vector<vector<double>>& house, vector<double> v, double b)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j) house[i][j] = 1.0 - b * v[i] * v[j];
			else house[i][j] = 0 - b * v[i] * v[j];
		}
	}
	/*cout << endl << endl;
	Matrix_cout(house, n, n);
	cout << endl << endl;*/
}

//局部矩阵提取+局部向量b；
void part_make(vector<vector<double>> matrix, int m,int n,vector<vector<double>>& part, vector<double> b, int k)
{
	for (int i = k; i < m; i++)
	{
		for (int j = k; j < n; j++)
			part[i - k][j - k] = matrix[i][j];
		part[i - k][n - k] = b[i];
	}
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

//Householder方法计算QR分解，matrix下方记录Householder生成向量，d记录每个Householder生成的系数值；
void QR(vector<vector<double>>& matrix,int m,int n,vector<double> &d)
{
	for (int j = 0; j < n; j++)
	{
		if (j < m - 1)
		{
			double b = 0;
			vector<double> v(m - j);
			vector<double> x(m - j);
			for (int i = 0; i < m - j; i++)
				x[i] = matrix[j + i][j];
			Householder(x, v, b);
			vector<vector<double>> house(m - j, vector<double>(m - j, 0));
			Householder_make(m - j, house, v, b);
			vector<vector<double>> temp(m - j, vector<double>(n - j + 1, 0));
			part_make(matrix, m, n, temp, d, j);
			Matrix_times(house, temp, m - j, n - j + 1);
			/*for (int i = j; i < m; i++)//右下局部作Householder变换；
			{
				for (int k = j; k < n; k++)
				{
					double temp = 0;
					for (int s = 0; s < m - j; s++)
						temp += return_value(v, b, i - j, s) * matrix[j + s][k];
					matrix[i][k] = temp;
				}
			}*/
			//d[j] = b;
			for (int i = j; i < m; i++)
			{
				for (int k = j; k < n; k++)
					matrix[i][k] = temp[i - j][k - j];
				d[i] = temp[i - j][n - j];
			}
			for (int i = 1; i < m - j; i++)
				matrix[j + i][j] = v[i];
		}
	}
}

/*//辅助函数，QR分解后，储存了生成向量与系数，构建返回第j个Householder方阵作乘法时ik位置所需的元素；
double return_value_jik(vector<vector<double>> matrix, double b, int j,int i, int k)
{
	if (i == k)
	{
		if (i == 0) return (1 - b);
		else return (1.0 - b * matrix[j][j + i] * matrix[j][j + k]);
	}
	else
	{
		if (i == 0) return (-b * 1.0 * matrix[j][j + k]);
		else if (k == 0) return (-b * matrix[j][j + i] * 1.0);
		else return (0 - b * matrix[j][j + i] * matrix[j][j + k]);
	}
}

//QR分解后的Q转置乘以向量b；
void Qt_times_b(vector<vector<double>> matrix, int m, int n, vector<double> d,vector<double>& b)
{
	for (int j = 0; j < n; j++)        //第j个Householder方阵；
	{
		for (int i = j; i < m; i++)    //b从第j个值开始变化，故乘法从第j行开始；
		{
			double temp = 0;
			for (int k = j; k < m; k++)//第i行乘以b得到新b[i]；
				temp += return_value_jik(matrix, d[j], j, i - j, k - j) * b[k];
			b[i] = temp;
		}
	}
}*/

//应用Householder阵作正交变换法求解LS问题，返回最小解x；
void LS(vector<vector<double>>& matrix, int m, int n, vector<double> &b)
{
	//vector<double> d(m);//向量d储存Householder变换的生成系数；
	QR(matrix, m, n, b);
	/*Qt_times_b(matrix, m, n, d, b);
	for (int i = 0; i < n; i++)
		x[i] = b[i];*/
	upper_triangular(n, matrix, b);
}

//2—范数的平方；
double norm_2(vector<vector<double>> matrix, int m, int n, vector<double> b, vector<double> x)
{
	double max = 0, temp;
	for (int i = 0; i < m; i++)
	{
		temp = 0;
		for (int j = 0; j < n; j++)
			temp += matrix[i][j] * x[j];
		max += pow(b[i] - temp, 2);
	}
	return max;
}