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
			cout << setw(9) << matrix[i][j] << " ";
		cout << endl;
	}
}

//矩阵乘法,并储存在第一个矩阵里面；
void Matrix_times(vector<vector<double>> &house, vector<vector<double>> matrix)
{
	int n = size(house);
	vector<vector<double>> temp(n, vector<double>(n, 0));
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				temp[i][j] += house[i][k] * matrix[k][j];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			house[i][j] = temp[i][j];
}

//矩阵转置；
vector<vector<double>> Matrix_transpose(vector<vector<double>> A)
{
	int dim = size(A);
	vector<vector<double>> T(dim, vector<double>(dim, 0));
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			T[i][j] = A[j][i];
	return T;
}

//返回矩阵的模最大元素；
double max_element_of_module(vector<vector<double>> A)
{
	int m = size(A);
	int n = size(A[0]);
	double max = fabs(A[0][0]);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			if (fabs(A[i][j]) > max)
				max = fabs(A[i][j]);
	return max;
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

//算法7.6.1辅助函数，Householder左作用；
void householder_left_operate(int m, int n, vector<vector<double>>& A, int k, vector<double> v, double b)
{
	int vlength = size(v);
	int ulength = n - k;
	vector<double> u(ulength, 0);
	for (int j = 0; j < ulength; j++)
		for (int i = 0; i < vlength; i++)
			u[j] += b * v[i] * A[k + i][k + j];
	for (int i = 0; i < vlength; i++)
		for (int j = 0; j < ulength; j++)
			A[k + i][k + j] -= v[i] * u[j];
}

//算法7.6.1辅助函数，Householder右作用；
void household_right_operate(int m, int n, vector<vector<double>>& A, int k, vector<double> v, double b)
{
	int vlength = size(v);
	int ulength = m - k;
	vector<double> u(ulength, 0);
	for (int i = 0; i < ulength; i++)
		for (int j = 0; j < vlength; j++)
			u[i] += A[k + i][k + 1 + j] * b * v[j];
	for (int i = 0; i < ulength; i++)
		for (int j = 0; j < vlength; j++)
			A[i + k][k + 1 + j] -= u[i] * v[j];
}

//算法7.6.1 Householder变换法二对角化，正交阵分别为U，V；
void bidiagonal_matrix(int m, int n, vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V)
{
	for (int k = 0; k < n; k++)
	{
		vector<double> x(m - k, 0), v(m - k, 0);
		double b = 0;
		for (int i = k; i < m; i++) x[i - k] = A[i][k];
		Householder(x, v, b);
		householder_left_operate(m, n, A, k, v, b);
		householder_left_operate(m, m, U, k, v, b);
		if (k < n - 2)
		{
			vector<double> y(n - k - 1, 0), u(n - k - 1, 0);
			b = 0;
			for (int j = k + 1; j < n; j++) y[j-k-1] = A[k][j];
			Householder(y, u, b);
			household_right_operate(m, n, A, k, u, b);
			household_right_operate(n, n, V, k, u, b);
		}
	}
	//清理相关元素；
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			if (fabs(A[i][j]) < 1e-10)
				A[i][j] = 0;
}

//算法7.6.3 Step3 收敛性判定；
void convergence_determination(int n, vector<vector<double>>& A, int& p, int & q)
{
	double epsilon = 1e-10;
	//求B的无穷范数；
	double infinity_norm = fabs(A[n - 1][n - 1]);
	for (int i = 0; i < n - 1; i++)
		if ((fabs(A[i][i]) + fabs(A[i][i + 1])) > infinity_norm)
			infinity_norm = fabs(A[i][i]) + fabs(A[i][i + 1]);
	//(i)；
	for (int i = 0; i < n - 1; i++)
		if (fabs(A[i][i + 1]) <= epsilon * (fabs(A[i][i]) + fabs(A[i + 1][i + 1])))
			A[i][i + 1] = 0;
	//(ii)；
	for (int i = 0; i < n; i++)
		if (fabs(A[i][i]) <= epsilon * infinity_norm)
			A[i][i] = 0;
	//(iii)确定最大非负整数p和最小非负整数q；
	p = 0, q = 0;
	int i = n - 2;
	while ((i >= 0) && (A[i][i + 1] == 0)) i--, q++;
	if (q == n - 1) q = n, p = 0;
	else
	{
		while ((i >= 0) && (A[i][i + 1] != 0)) i--;
		p = i + 1;
	}
	/*int qflag = n - 2;
	while ((qflag >= 0) && (A[qflag][qflag + 1] == 0)) qflag--, q++;
	if (q < n)
	{
		while ((qflag >= 0) && (A[qflag][qflag] != 0)) qflag--;
		p = qflag + 1;
	}*/
}

//Givens变换；
void Givens_ab(double a, double b, double& c, double& s)
{
	double t;
	if (b == 0) c = 1, s = 0;
	else
		if (abs(b) > abs(a)) t = a / b, s = 1 / sqrt(1 + t * t), c = s * t;
		else t = b / a, c = 1 / sqrt(1 + t * t), s = c * t;
}

//算法7.6.3 Step4(i) 判定B22对角元是否有0元素，其中fla用于是否有，有1、无0；
void svd_iteration_determination(int m, int n, vector<vector<double>>& A, vector<vector<double>>&  U, int p, int q, int& flag)
{
	flag = 0;
	for (int i = p; i < n; i++)
	{
		if (A[i][i] == 0) 
		{
			flag = 1;
			double c = 0, s = 0;
			int k = 1;
			while (i + k <= n - 1)
			{
				Givens_ab(A[i][i + k], A[i + k][i + k], c, s);
				//记录A变化；
				double temp = A[i][i + k];
				A[i][i + k] = c * A[i][i + k] + s * A[i + k][i + k];//实质为0；
				//A[i][i + k] = 0;
				A[i + k][i + k] = -s * temp + c * A[i + k][i + k];
				if (i + k < n - 1)
				{
					A[i][i + k + 1] = s * A[i + k][i + k + 1];
					A[i + k][i + k + 1] = c * A[i + k][i + k + 1];
				}
				//记录U变化；
				vector<double> vi(m, 0), vik(m, 0);//记录两行变化；
				for (int j = 0; j < m; j++) 
					vi[j] = c * U[i][j] + s * U[i + k][j], vik[j] = -s * U[i][j] + c * U[i + k][j];
				for (int j = 0; j < m; j++)
					U[i][j] = vi[j], U[i + k][j] = vik[j];
				//i+k累加至n-1；
				k++;
			}
		}
	}
}

//算法7.6.2 带Wilkinson位移的SVD迭代；
void Wilkinson_svd_iteration(vector<vector<double>>& B, vector<vector<double>>& U, vector<vector<double>>& V, int b)
{
	int p = size(U), q = size(V), n = size(B);
	double t = 0;
	if (n >= 3) t = B[n - 3][n - 2];
	double alpha = pow(B[n - 1][n - 1], 2) + pow(B[n - 2][n - 1], 2);
	double delta = (pow(B[n - 2][n - 2], 2) + pow(t, 2) - alpha) / 2;
	double beta = B[n - 2][n - 2] * B[n - 2][n - 1];
	double sign = 1.0;
	if (delta < 0) sign = -1.0;
	double miu = alpha - pow(beta, 2) / (delta + sign * sqrt(pow(delta, 2) + pow(beta, 2)));
	double y = pow(B[0][0], 2) - miu;
	double z = B[0][0] * B[0][1];
	for (int k = 1; k < n; k++)
	{
		double c = 0, s = 0;
		c = y / (sqrt(y * y + z * z)), s = -z / (sqrt(y * y + z * z));
		if (k > 1) B[k - 2][k - 1] = c * y - s * z;
		y = B[k - 1][k - 1] * c - B[k - 1][k] * s;
		z = -s * B[k][k];
		B[k - 1][k] = B[k - 1][k - 1] * s + B[k - 1][k] * c;
		B[k][k] = c * B[k][k];
		//记录Q于V；
		double u = 0, v = 0;
		for (int i = 0; i < q; i++)
		{
			u = V[i][b + k - 1] * c - V[i][b + k] * s;
			v = V[i][b + k - 1] * s + V[i][b + k] * c;
			V[i][b + k - 1] = u, V[i][b + k] = v;
		}

		c = y / (sqrt(y * y + z * z)), s = -z / (sqrt(y * y + z * z));
		B[k - 1][k - 1] = c * y - s * z;
		if (k < n - 1)
		{
			y = c * B[k - 1][k] - s * B[k][k];
			z = -s * B[k][k + 1];
			B[k][k] = s * B[k - 1][k] + c * B[k][k];
			B[k][k + 1] = c * B[k][k + 1];
		}
		else B[k - 1][k] = c * B[k - 1][k] - s * B[k][k], B[k][k] = s * B[k - 1][k] + c * B[k][k];
		//记录P于U；
		for (int i = 0; i < p; i++)
		{
			u = U[b + k - 1][i] * c - U[b + k][i] * s;
			v = U[b + k - 1][i] * s + U[b + k][i] * c;

			//u = c * U[b + k - 1][i] + s * U[b + k][i];
			//v = -s * U[b + k - 1][i] + c * U[b + k][i];
			U[b + k - 1][i] = u, U[b + k][i] = v;
		}
	}
}

//算法7.6.3SVD分解；
void SVD_decomposition(int m, int n, vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V)
{
	//Step2；
	bidiagonal_matrix(m, n, A, U, V);
	//Step3；
	int p = 0, q = 0, flag = 0;
	convergence_determination(n, A, p, q);
	while (q != n)
	{
		//Step4(i)；
		svd_iteration_determination(m, n, A, U, p, q, flag);
		//cout << endl << flag << endl;
		while (flag == 1)
		{
			convergence_determination(n, A, p, q);
			svd_iteration_determination(m, n, A, U, p, q, flag);
		}
		if (q != n)//判定终止条件；
		{
			//提取B22；
			int B22_dim = n - p - q;
			vector<vector<double>> B22(B22_dim, vector<double>(B22_dim, 0));
			for (int i = 0; i < B22_dim; i++)
				for (int j = 0; j < B22_dim; j++)
					B22[i][j] = A[i + p][j + p];
			//Step4(ii) 带Wilkinson位移的SVD迭代；
			Wilkinson_svd_iteration(B22, U, V, p);
			//赋值回去；
			for (int i = 0; i < B22_dim; i++)
				for (int j = 0; j < B22_dim; j++)
					A[i + p][j + p] = B22[i][j];
			//转Step3；
			convergence_determination(n, A, p, q);
		}
	}
}

//返回奇异值分解乘机；
vector<vector<double>> UBV(vector<vector<double>> U, vector<vector<double>> A, vector<vector<double>> V)
{
	int m = size(U);
	int n = size(V);
	vector<vector<double>> B(m, vector<double>(n, 0));
	vector<vector<double>> temp(m, vector<double>(n, 0));
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < m; k++)
				B[i][j] += U[i][k] * A[k][j];
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				temp[i][j] += B[i][k] * V[k][j];
	return temp;
}