#include <iostream>
#include<iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <ctime>

using namespace std;

//矩阵输出；
void Matrix_cout(vector<vector<double>> matrix, int m, int n)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			cout << setw(9) << matrix[i][j] << "    ";
		cout << endl;
	}
}

//矩阵乘法house*matrix,并储存在第三个矩阵里面；
void Matrix_times(vector<vector<double>> house, vector<vector<double>> matrix, vector<vector<double>>& A)
{
	int dim = size(matrix);
	vector<vector<double>> temp(dim, vector<double>(dim, 0));
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			for (int k = 0; k < dim; k++)
				temp[i][j] += house[i][k] * matrix[k][j];
		}
	}
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			A[i][j] = temp[i][j];
}

//生成多项式（向量储存）友方阵；
void Companion_matrix(vector<vector<double>>& matrix, vector<double> poly)
{
	int dim = size(poly);
	int n = dim - 1;
	for (int i = 0; i < n; i++)
		matrix[i][n - 1] = 0 - poly[i];
	for (int i = 1; i < n; i++)
		matrix[i][i - 1] = 1.0;
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

//向量的模最大分量；
double magnitude_max_component_of_vector(vector<double> x)
{
	int dim = size(x);
	double norm = fabs(x[0]);
	int temp = 0;
	for (int i = 1; i < dim; i++)
	{
		if (fabs(x[i]) > norm)
		{
			norm = fabs(x[i]);
			temp = i;
		}
	}
	return x[temp];
}

//向量之间作减法；
void Vector_subtraction(vector<double> x, vector<double> y, vector<double>& dx)
{
	int dim = 0;
	dim = size(x);
	for (int i = 0; i < dim; i++)
		dx[i] = x[i] - y[i];
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

//辅助函数，由yk，muk生成uk；
void uk(vector<double> yk, double muk, vector<double>& uk)
{
	int dim = size(yk);
	for (int i = 0; i < dim; i++)
		uk[i] = yk[i] / muk;
}

//幂法求解方阵模最大特征值；
void Power_method(int dim, vector<vector<double>> A, double error , double& value)
{
	vector<double> u0(dim, 1);//初始向量；
	vector<double> y(dim, 0), u(dim, 1), du(dim, 1);
	double mu = 0;
	do {
		Matrix_times_vector(A, u0, y);//求解向量yk
		mu = magnitude_max_component_of_vector(y);//求yk模最大分量；
		uk(y, mu, u);
		Vector_subtraction(u, u0, du);
		u0 = u;
	} while (vector_infinity_norm(du) > error);
	//value = mu;//利用迭代得出的mu给出模最大特征值，普遍适用；
	value = A[0][dim - 1] * u[dim - 1] / u[0];//利用特征向量与友方阵特性，以及特征值定义，利用第一个分量求解特征值，只适用友方阵情形；
}

//幂法作用于多项式友方阵求解模最大根；
double max_modulus_root_polynomial(vector<double> poly, double error)
{
	int dim = size(poly);
	int n = dim - 1;//多项式阶数，矩阵维数；
	vector<vector<double>> matrix(n, vector<double>(n, 0));
	//生成友方阵；
	for (int i = 0; i < n; i++)
		matrix[i][n - 1] = 0 - poly[dim - 1 - i];
	//for (int i = 0; i < n; i++)//老师上课讲的友方阵，实际只是求出了一个特征值，但不是模最大；
	//	matrix[0][i] = 0 - poly[i + 1];
	for (int i = 1; i < n; i++)
		matrix[i][i - 1] = 1.0;
	//Matrix_cout(matrix, n, n);
	double value = 0;
	Power_method(n, matrix, error, value);
	return value;
}

//计算Householder变换，使得向量变为e1，返回生成向量v，系数b；
void Householder(vector<double> x, vector<double>& v, double& b)
{
	int dim = size(x);
	double max = 0, m = 0, a = 0;
	max = vector_infinity_norm(x);
	x[0] = x[0] / max;
	for (int i = 1; i < dim; i++)
	{
		v[i] = x[i] / max;
		m += pow(v[i], 2);
	}
	if (m == 0) b = 0;
	else
	{
		a = sqrt(x[0] * x[0] + m);
		if (x[0] <= 0) v[0] = x[0] - a;
		else v[0] = -m / (x[0] + a);
		b = 2 * pow(v[0], 2) / (m + pow(v[0], 2));
		for (int i = 1; i < dim; i++)
			v[i] = v[i] / v[0];
		v[0] = 1.0;
	}
}

//Householder方阵生成；
void Householder_make(int n, vector<vector<double>>& house, vector<double> v, double b)
{
	int dim = size(house);
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			if (i > n && j > n) 
			{
				if (i == j) house[i][j] = 1.0 - b * v[i - n - 1] * v[j - n - 1];
				else house[i][j] = 0 - b * v[i - n - 1] * v[j - n - 1];
			}
			else
			{
				if (i == j) house[i][j] = 1;
				else house[i][j] = 0;
			}
		}
	}
}

//6.4.1辅助函数，Householder阵左右作用；
void householder_opertion(vector<vector<double>>& A, int n, vector<double> v,double b)
{
	int dim = size(A);
	int vl = size(v);
	vector<vector<double>> temp(vl, vector<double>(vl, 0));
	for (int i = 0; i < vl; i++)
		for (int j = 0; j < vl; j++)
			for (int k = 0; k < vl; k++)
				temp[i][j] += b * v[i] * v[j] * A[n + 1 + k][n + 1 + j];
	for (int i = 0; i < vl; i++)
		for (int j = 0; j < vl; j++)
			A[n + 1 + i][n + 1 + j] -= temp[i][j];
	for (int i = 2; i < dim; i++)
		for (int j = 0; j < i - 1; j++)
			A[i][j] = 0;

	vector<vector<double>> temp0(vl, vector<double>(vl, 0));
	for (int i = 0; i < vl; i++)
		for (int j = 0; j < vl; j++)
			for (int k = 0; k < vl; k++)
				temp0[i][j] += A[n + 1 + i][n + 1 + k] * b * v[k] * v[j];
	for (int i = 0; i < vl; i++)
		for (int j = 0; j < vl; j++)
			A[n + 1 + i][n + 1 + j] -= temp0[i][j];
	for (int i = 2; i < dim; i++)
		for (int j = 0; j < i - 1; j++)
			A[i][j] = 0;
}

//算法6.4.1 用Householder方法计算矩阵的上Hessenberg分解，分解储存在A里，并将变换矩阵储存在Q里；
void Hessenberg_decomposition(vector<vector<double>>& A)
{
	int dim = size(A);
	for (int j = 0; j < dim - 2; j++)
	{
		double b = 0;
		int tempk = dim - j - 1;
		vector<double> v(tempk, 0);
		vector<double> x(tempk, 0);  //提取i列对角元以下用于Householder变换；
		for (int i = 0; i < tempk; i++)
			x[i] = A[i + j + 1][j];
		Householder(x, v, b);        //求Householder变换；
		vector<vector<double>> house(dim, vector<double>(dim, 0));
		Householder_make(j, house, v, b);//生成Householder变换阵；
		Matrix_times(house, A, A);//Householder变换左乘；
		for (int i = j + 2; i < dim; i++) A[i][j] = 0;
		Matrix_times(A, house, A);//Householder变换右乘；
		
	}
	//清理下三角阵；
	for (int i = 2; i < dim; i++)
		for (int j = 0; j < i - 1; j++)
			A[i][j] = 0;
}

//清理Hesssenberg阵中小元素置为0；
void clean_small_elements(int dim, vector<vector<double>>& H, double u)
{
	for (int i = 1; i < dim; i++)
		if (fabs(H[i][i - 1]) <= (fabs(H[i][i]) + fabs(H[i - 1][i - 1])) * u) H[i][i - 1] = 0;
}

//确定隐式QR中提到的m,l参数，注意m为已有的m阶拟上三角阵，默认H是上Hessenberg阵
void confirm_m_l(vector<vector<double>> H, int& m, int& l)
{
	int n = H.size(), t = 1;
	while ((m < n - 2) && (t == 1))
	{
		if (H[n - m - 1][n - m - 2] == 0) m++;
		else if ((H[n - m - 2][n - m - 3] == 0) && ((H[n - m - 2][n - m - 2] - H[n - m - 1][n - m - 1]) * (H[n - m - 2][n - m - 2] - H[n - m - 1][n - m - 1]) + 4 * H[n - m - 2][n - m - 1] * H[n - m - 1][n - m - 2]) < 0)
			//判定2阶情形
		{
			m += 2, t += 2;
		}
		else t = 0;
	}
	if (m == n - 1) m = n, l = 0;
	else if (m == n - 2)
	{
		if (((H[0][0] - H[1][1]) * (H[0][0] - H[1][1]) + 4 * H[0][1] * H[1][0] < 0)) m = n, l = 0;
		else l = 2;
	}

	if (m != n) for (l = 2; (H[n - m - l + 1][n - m - l] != 0) && (l < (n - m)); l++);
	l = n - m - l;
}

//双重步位移的QR迭代&默认P为0阵
void Francis(vector<vector<double>>& H, vector<vector<double>>& P)
{
	int n = H.size(), q, r, i;
	double s, t, x, y, z, a, b, c, d;
	vector<double> v(3), w(3);
	for (i = 0; i < n; i++) P[i][i] = 1;
	if (n >= 3)
	{
		s = H[n - 2][n - 2] + H[n - 1][n - 1];
		t = H[n - 2][n - 2] * H[n - 1][n - 1] - H[n - 2][n - 1] * H[n - 1][n - 2];
		x = H[0][0] * H[0][0] + H[0][1] * H[1][0] - s * H[0][0] + t;
		y = H[1][0] * (H[0][0] + H[1][1] - s);
		z = H[1][0] * H[2][1];
		for (int k = 0; k < n - 2; k++)
		{
			a = sqrt(x * x + y * y + z * z);
			v[0] = x - a, v[1] = y, v[2] = z;
			b = a * a - a * x;                       //简化Householder变换中的beta表达式；
			if (b != 0)
			{
				b = 1 / b;
				q = (1 > k) ? 1 : k;
				w[0] = b * v[0], w[1] = b * v[1], w[2] = b * v[2];
				for (i = q - 1; i < n; i++)
				{
					c = v[0] * H[k][i] + v[1] * H[k + 1][i] + v[2] * H[k + 2][i];
					H[k][i] -= w[0] * c, H[k + 1][i] -= w[1] * c, H[k + 2][i] -= w[2] * c;
				}
				r = (k + 4 > n) ? n : k + 4;
				for (i = 0; i < r; i++)
				{
					c = H[i][k] * v[0] + H[i][k + 1] * v[1] + H[i][k + 2] * v[2];
					H[i][k] -= c * w[0], H[i][k + 1] -= c * w[1], H[i][k + 2] -= c * w[2];
					d = P[i][k] * v[0] + P[i][k + 1] * v[1] + P[i][k + 2] * v[2];
					P[i][k] -= d * w[0], P[i][k + 1] -= d * w[1], P[i][k + 2] -= d * w[2];
				}
				x = H[k + 1][k], y = H[k + 2][k];
				if (k < n - 3) z = H[k + 3][k];
			}
		}
		a = sqrt(x * x + y * y);
		v[0] = x - a, v[1] = y;
		b = a * a - a * x;
		if (b != 0)
		{
			b = 1 / b;
			w[0] = b * v[0], w[1] = b * v[1];
			for (i = n - 3; i < n; i++)
			{
				c = v[0] * H[n - 2][i] + v[1] * H[n - 1][i];
				H[n - 2][i] -= w[0] * c, H[n - 1][i] -= w[1] * c;
			}
			for (i = 0; i < n; i++)
			{
				c = H[i][n - 2] * v[0] + H[i][n - 1] * v[1];
				H[i][n - 2] -= w[0] * c, H[i][n - 1] -= w[1] * c;
				d = P[i][n - 2] * v[0] + P[i][n - 1] * v[1];
				P[i][n - 2] -= w[0] * d, P[i][n - 1] -= w[1] * d;
			}
		}
	}
	else if (n == 2)
		//能通过Judgegment函数进来的方阵必然满足s非负，所以P不可能为单位阵
	{
		a = H[0][0], b = H[0][1], c = H[1][0], d = H[1][1];
		s = (a - d) * (a - d) + 4 * b * c;
		if (s >= 0)
		{
			H[0][1] = b - c, H[1][0] = 0, H[0][0] = (a + d - sqrt(s)) / 2, H[1][1] = (a + d + sqrt(s)) / 2;
		}
		if (b == 0) P[0][0] = 0, P[0][1] = 1, P[1][0] = -1, P[1][1] = 0;
		else t = (a - H[0][0]) / b, P[0][0] = P[1][1] = 1 / sqrt(1 + t * t), P[0][1] = P[0][0] * t, P[1][0] = -P[0][1];
	}
}

//6.4.3辅助函数，提取矩阵；
void part_make(vector<vector<double>> H, vector<vector<double>>& part, int starti, int endi, int startj, int endj)
{
	for (int i = starti; i <= endi; i++)
		for (int j = startj; j <= endj; j++)
			part[i - starti][j - startj] = H[i][j];
}

//6.4.3辅助函数，H12和H23乘P的操作；
void submatrix_operations(vector<vector<double>>& H, vector<vector<double>> P, int m, int l)
{
	int n = H.size();
	vector<vector<double>> Q(l, vector<double>(n - m - l, 0)), R(n - m - l, vector<double>(m, 0));
	//H12的操作
	for (int i = 0; i < l; i++)
		for (int j = 0; j < n - m - l; j++)
			for (int k = 0; k < n - m - l; k++)
				Q[i][j] += H[i][k + l] * P[k][j];
	for (int i = 0; i < l; i++)
		for (int j = 0; j < n - m - l; j++)
			H[i][j + l] = Q[i][j];
	//H23的操作
	for (int i = 0; i < n - m - l; i++)
		for (int j = 0; j < m; j++)
			for (int k = 0; k < n - m - l; k++)
				R[i][j] += P[k][i] * H[l + k][n - m + j];
	for (int i = 0; i < n - m - l; i++)
		for (int j = 0; j < m; j++)
			H[l + i][n - m + j] = R[i][j];
}

//算法6.4.3 计算实矩阵的实Schur分解，即隐式QR算法；
void Real_Schur_decomposition(vector<vector<double>>& A, double u)
{
	int dim = size(A);
	int m = 0, l = 0;
	Hessenberg_decomposition(A);        //Step2；
	clean_small_elements(dim, A, u);    //Step3(1)；
	confirm_m_l(A, m, l);
	while (m != dim)
	{
		//Step4；
		vector<vector<double>> H22(dim - m - l, vector<double>(dim - m - l, 0)), P(dim - m - l, vector<double>(dim - m - l, 0));
		part_make(A, H22, l, dim - m - 1, l, dim - m - 1);//取出H22；
		Francis(H22, P);
		for (int i = l; i < dim - m; i++)                 //赋值回去；
			for (int j = l; j < dim - m; j++)
				A[i][j] = H22[i - l][j - l];

		//Step5；
		submatrix_operations(A, P, m, l);
		for (int i = 1; i < dim; i++) 
			if (abs(A[i][i - 1]) < 0.00001) A[i][i - 1] = 0;

		//清理下三角阵；
		for (int i = 2; i < dim; i++)
			for (int j = 0; j < i - 1; j++)
				A[i][j] = 0;
		
		clean_small_elements(dim, A, u);//Step3(1)；
		confirm_m_l(A, m, l);           //Step3(2)；
	}
	for (int i = 1; i < dim; i++) if (abs(A[i][i - 1]) < 0.0000001) A[i][i - 1] = 0;
	
	/*用do while 写法，无本质差异,但更符合6.4.3的顺序；
	do {
		clean_small_elements(dim, A, u);//Step3(1)；
		int m = 0, l = 0;
		confirm_m_l(A, m, l);           //Step3(2)；
		if (m >= dim) break;            //Step3(3)；

		//Step4；
		vector<vector<double>> H22(dim - m - l, vector<double>(dim - m - l, 0)), P(dim - m - l, vector<double>(dim - m - l, 0));
		part_make(A, H22, l, dim - m - 1, l, dim - m - 1);//取出H22；
		Francis(H22, P);
		for (int i = l; i < dim - m; i++)                  //赋值回去；
			for (int j = l; j < dim - m; j++)
				A[i][j] = H22[i - l][j - l];

		//Step5；
		submatrix_operations(A, P, m, l);
		for (int i = 1; i < dim; i++) if (abs(A[i][i - 1]) < 0.00001) A[i][i - 1] = 0;

		//清理下三角阵；
		for (int i = 2; i < dim; i++)
			for (int j = 0; j < i - 1; j++)
				A[i][j] = 0;
	} while (1);*/
}

//给定方阵的Schur分解，求解方阵的特征值并输出；
void Eigenvalue(vector<vector<double>> A)
{
	int dim = size(A);
	int i = 0;
	while (i < dim - 1)
	{
		if (A[i + 1][i] == 0)
			cout << setw(9) << A[i][i] << "  ", i++;
		else 
		{
			double a = A[i][i], b = A[i][i + 1], c = A[i + 1][i], d = A[i + 1][i + 1];
			cout << setw(9) << (a + d) / 2 << "±" << sqrt(-4 * b * c - (a - d) * (a - d)) / 2 << "i" << "  ", i += 2;
		}
	}
	if (i == dim - 1)
		cout << setw(9) << A[dim - 1][dim - 1] << "  ";
	cout << endl;
}