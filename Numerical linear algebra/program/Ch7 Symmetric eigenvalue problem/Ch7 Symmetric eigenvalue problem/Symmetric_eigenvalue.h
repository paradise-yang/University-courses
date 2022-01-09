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

//非对角范数值；
double off_diagonal_norm(vector<vector<double>> A)
{
	int dim = size(A);
	double temp = 0;
	for (int i = 1; i < dim; i++)
		temp += 2 * A[i][i - 1] * A[i][i - 1];
	return sqrt(temp);
}

//确定Givens变换阵的值；
void Givens(vector<vector<double>> A, int p, int q, double& c, double& s)
{
	if (A[p][q] == 0) c = 1, s = 0;
	else 
	{
		int dim = size(A);
		double tao = (A[q][q] - A[p][p]) / (2 * A[p][q]);
		int sign = 1;//其中tao=0时候，sign=1,为了兼容；
		if (tao < 0) sign = -1.0;
		double t = sign / (fabs(tao) + sqrt(1 + tao * tao));
		c = 1 / sqrt(1 + t * t);
		s = t * c;
	}
}

//一次Jacobi扫描；
void one_Jacobi_scan(vector<vector<double>>& A, double deta)
{
	int dim = size(A);
	for (int p = 0; p < dim - 1; p++)
	{
		for (int q = p + 1; q < dim; q++)
		{
			if (fabs(A[p][q]) > deta)
			{
				
				double c = 0, s = 0;
				Givens(A, p, q, c, s);

				vector<double> vp(dim, 0), vq(dim, 0);//记录两行两列变化；
				for (int i = 0; i < dim; i++) vp[i] = c * A[i][p] - s * A[i][q];
				vp[p] = c * c * A[p][p] - 2 * s * c * A[p][q] + s * s * A[q][q];
				for (int i = 0; i < dim; i++) vq[i] = s * A[i][p] + c * A[i][q];
				vq[q] = s * s * A[p][p] + 2 * s * c * A[p][q] + c * c * A[q][q];
				vp[q] = vq[p] = (c * c - s * s) * A[p][q] + s * c * (A[p][p] - A[q][q]);

				for (int i = 0; i < dim; i++) A[i][p] = A[p][i] = vp[i], A[i][q] = A[q][i] = vq[i];
				//A[p][q] = A[q][p] = 0;
			}
		}
	}
}

//过关Jacobi方法求解矩阵特征值；
void Pass_Jacobi_method(vector<vector<double>>& A, double u)
{
	int dim = size(A);//sigma也取dim；
	double EA = off_diagonal_norm(A);//非对角范数；
	//cout << "第 0 次扫描后，当前非对角范数为 ：" << EA << endl;
	double deta = EA;
	while(deta > u)
	{
		deta = deta / dim;
		one_Jacobi_scan(A, deta);
		//cout << "第 " << i << " 次扫描后，当前非对角范数为 ：" << off_diagonal_norm(A) << endl;
	}
	cout << dim << " 阶矩阵特征值为：" << endl;
	for (int i = 0; i < dim; i++) cout << setw(9) << A[i][i] << " ";
	cout << endl;
}

//算法7.4.1计算变号数；
int variable_number(vector<vector<double>> A, double miu)
{
	int dim = size(A);
	double u = 1e-10;
	vector<double> x(dim, 0), y(dim, 0);
	for (int i = 0; i < dim; i++) x[i] = A[i][i];
	for (int i = 1; i < dim; i++) y[i] = A[i][i - 1];
	int s = 0;
	double q = x[1] - miu;
	for (int k = 0; k < dim; k++)
	{
		if (q < 0) s++;
		if (k < dim - 1)
		{
			if (q == 0) q = fabs(y[k + 1]) * u;
			q = x[k + 1] - miu - y[k + 1] * y[k + 1] / q;
		}
	}
	return s;
}

//二分法求指定特征值；
double bisection_method(vector<vector<double>> A, int k)
{
	double length = 1e-10;
	int dim = size(A);
	double max = 4;
	/*double max = fabs(A[0][0]) + fabs(A[0][1]);
	for (int i = 1; i < dim; i++)
	{
		double temp = fabs(A[i][i]) + 2 * fabs(A[i][i - 1]);
		if (temp > max) max = temp;
	}//确定无穷范数；
	*/
	double low = -max, up = max;
	while ((up - low) > length)
	{
		double middle = (low + up) / 2;
		if (variable_number(A, middle) >= k) up = middle;
		else low = middle;//, k = k - s;
	}
	return (low + up) / 2;
}