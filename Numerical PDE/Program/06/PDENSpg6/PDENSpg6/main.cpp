#include<iostream>
#include<iomanip>
#include<vector>
#include<cmath>
#include<fstream>
#include<string>

using namespace std;

#define pi 3.1415926535897932384626433832795028841971 //圆周率\pi值；
#define miu 0.1 //网格比；

//精确解函数；
double u_exact(double x, double t)
{
	if (x >= -pi && x <= 0)
		return 0.5*exp(-4 * t)*sin(x);
	else
		return exp(-4 * t)*sin(2 * x);
}
//初值函数；
double a(double x)
{
	if (x < 0)
		return 4.0;
	else
		return 1.0;
}
//FTBS格式；
void FTCS_mean(int J, double lamuda, double T, vector<vector<double>> &u)
{
	double deltax = 2.0 * pi / J;              //空间离散程度；
	double deltat = lamuda * deltax * deltax;//求解空间离散程度；
	int N = floor(T / deltat);               //需要前进的步数；
	double last_deltat = T - N * deltat;
	//初值；
	for (int j = 0; j <= J; j++)
		u[0][j] = u_exact(j * deltax - pi, 0);
	//FTBS格式；
	for (int i = 1; i <= N; i++)
	{
		u[i % 2][0] = 0;//左边界处理；
		for (int j = 1; j < J; j++)
		{
			double aleft = (a((j - 1) * deltax - pi) + a(j * deltax - pi)) / 2;
			double aright= (a((j + 1) * deltax - pi) + a(j * deltax - pi)) / 2;
			u[i % 2][j] = lamuda * aleft*u[(i + 1) % 2][j - 1] + (1 - lamuda * aright - lamuda * aleft)*u[(i + 1) % 2][j] + lamuda * aright*u[(i + 1) % 2][j + 1];
		}
		u[i % 2][J] = 0;//右边界处理；
	}
	//最后一个短步长处理；
	double lamuda_last = last_deltat / deltax / deltax;
	int i = N + 1;
	u[i % 2][0] = 0;//左边界处理；
	for (int j = 1; j < J; j++)
	{
		double aleft = (a((j - 1) * deltax - pi) + a(j * deltax - pi)) / 2;
		double aright = (a((j + 1) * deltax - pi) + a(j * deltax - pi)) / 2;
		u[i % 2][j] = lamuda * aleft*u[(i + 1) % 2][j - 1] + (1 - lamuda * aright - lamuda * aleft)*u[(i + 1) % 2][j] + lamuda * aright*u[(i + 1) % 2][j + 1];
	}
	u[i % 2][J] = 0;//右边界处理；
}
//
double theta(double x, double xj_1, double xj, double deltax)
{
	if (x >= xj_1 && x <= xj)
		return 1.0 / 2;
		//return (x - xj) / deltax;
	else
		return (x - xj) / deltax;
		//return 1.0 / 2;
}
void FTCS_harmonic(int J, double lamuda, double T, vector<vector<double>> &u)
{
	double deltax = 2.0 * pi / J;              //空间离散程度；
	double deltat = lamuda * deltax * deltax;//求解空间离散程度；
	int N = floor(T / deltat);               //需要前进的步数；
	double last_deltat = T - N * deltat;
	//初值；
	for (int j = 0; j <= J; j++)
		u[0][j] = u_exact(j * deltax - pi, 0);
	//FTBS格式；
	for (int i = 1; i <= N; i++)
	{
		u[i % 2][0] = 0;//左边界处理；
		for (int j = 1; j < J; j++)
		{
			double aleft = 1.0 / (theta(0, (j - 1) * deltax - pi, j * deltax - pi, deltax) / a((j - 1) * deltax - pi) + (1.0 - theta(0, (j - 1) * deltax - pi, j * deltax - pi, deltax)) / a(j * deltax - pi));
				//(a((j - 1) * deltax - pi) + a(j * deltax - pi)) / 2;
			double aright = 1.0 / (theta(0, j * deltax - pi, (j + 1) * deltax - pi, deltax) / a(j * deltax - pi) + (1.0 - theta(0, j * deltax - pi, (j + 1) * deltax - pi, deltax)) / a((j + 1) * deltax - pi));
				//(a((j + 1) * deltax - pi) + a(j * deltax - pi)) / 2;
			u[i % 2][j] = lamuda * aleft*u[(i + 1) % 2][j - 1] + (1 - lamuda * aright - lamuda * aleft)*u[(i + 1) % 2][j] + lamuda * aright*u[(i + 1) % 2][j + 1];
		}
		u[i % 2][J] = 0;//右边界处理；
	}
	//最后一个短步长处理；
	double lamuda_last = last_deltat / deltax / deltax;
	int i = N + 1;
	u[i % 2][0] = 0;//左边界处理；
	for (int j = 1; j < J; j++)
	{
		double aleft = 1.0 / (theta(0, (j - 1) * deltax - pi, j * deltax - pi, deltax) / a((j - 1) * deltax - pi) + (1.0 - theta(0, (j - 1) * deltax - pi, j * deltax - pi, deltax)) / a(j * deltax - pi));
		double aright = 1.0 / (theta(0, j * deltax - pi, (j + 1) * deltax - pi, deltax) / a(j * deltax - pi) + (1.0 - theta(0, j * deltax - pi, (j + 1) * deltax - pi, deltax)) / a((j + 1) * deltax - pi));
		u[i % 2][j] = lamuda * aleft*u[(i + 1) % 2][j - 1] + (1 - lamuda * aright - lamuda * aleft)*u[(i + 1) % 2][j] + lamuda * aright*u[(i + 1) % 2][j + 1];
	}
	u[i % 2][J] = 0;//右边界处理；
}
//求解模最大误差；
double error_max(int J, double lamuda, double T, vector<vector<double>> u)
{
	double deltax = 2.0 * pi / J;              //空间离散程度；
	double deltat = lamuda * deltax * deltax;//求解空间离散程度；
	int N = floor(T / deltat) + 1;           //需要前进的步数；
	double max = 0;                 //记录误差；
	for (int j = 0; j <= J; j++)
		if (fabs(u_exact(j*deltax - pi, T) - u[N % 2][j]) > max)
			max = fabs(u_exact(j*deltax - pi, T) - u[N % 2][j]);
	return max;
}
//求解L2模误差；
double error_L2(int J, double lamuda, double T, vector<vector<double>> u)
{
	double deltax = 2.0 * pi / J;              //空间离散程度；
	double deltat = lamuda * deltax * deltax;//求解空间离散程度；
	int N = floor(T / deltat) + 1;           //需要前进的步数；
	double max = 0;                 //记录误差；
	for (int j = 0; j <= J; j++)
		max += pow(u_exact(j*deltax - pi, T) - u[N % 2][j], 2);
	return sqrt(max*deltax);
}
//矩阵输出；
void output_matrix(vector<vector<double>> matrix)
{
	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix[0].size(); j++)
			cout << matrix[i][j] << "\t";
		cout << endl;
	}
}
//输出数值解进入指定文本文件；
void output(int i, double deltax, double lamuda, double T, vector<vector<double>> u)
{
	double deltat = lamuda * deltax;//求解空间离散程度；
	int N = T / deltat;             //需要前进的步数；
	ofstream outfile;//输出三维数值误差；
	string temp = "E:\\study_materials\\PDEns\\06\\ns";
	string path = temp + to_string(i) + ".txt";
	outfile.open(path, ios::out);
	for (int j = 0; j < size(u[0]); j++)
		outfile << j * deltax << " " << u[N % 2][j] << endl;
}

int main()
{
	vector<vector<double>> error(5, vector<double>(5, 0));
	vector<vector<double>> error2(5, vector<double>(5, 0));
	for (int i = 1; i <= 5; i++)
	{
		int J = 10 * pow(2, i) + 1;                           //空间离散程度；
		double deltax = 2 * pi / J;                           
		double T = 1.0;                                       //待求时间T处数值解；
		vector<vector<double>> u(3, vector<double>(J + 1, 0));//记录求解向量；
		vector<vector<double>> u2(3, vector<double>(J + 1, 0));//记录求解向量；
		FTCS_mean(J, miu, T, u);
		FTCS_harmonic(J, miu, T, u2);
		error[i - 1][0] = J, error2[i - 1][0] = J;
		error[i - 1][1] = error_max(J, miu, T, u), error2[i - 1][1] = error_max(J, miu, T, u2);
		error[i - 1][3] = error_L2(J, miu, T, u), error2[i - 1][3] = error_L2(J, miu, T, u2);
	}
	for (int i = 1; i < 5; i++)
		error[i][2] = log(error[i - 1][1] / error[i][1]) / log((10 * pow(2, i + 1) + 1) / (10 * pow(2, i) + 1)),
		error[i][4] = log(error[i - 1][3] / error[i][3]) / log((10 * pow(2, i + 1) + 1) / (10 * pow(2, i) + 1)),
		error2[i][2] = log(error2[i - 1][1] / error2[i][1]) / log((10 * pow(2, i + 1) + 1) / (10 * pow(2, i) + 1)),
		error2[i][4] = log(error2[i - 1][3] / error2[i][3]) / log((10 * pow(2, i + 1) + 1) / (10 * pow(2, i) + 1));
	cout << "算术平均：" << endl;
	cout << "N\t模最大误差\t收敛阶\t模L2误差\t收敛阶" << endl;
	output_matrix(error);
	cout << endl << "调和平均：" << endl;
	cout << "N\t模最大误差\t收敛阶\t模L2误差\t收敛阶" << endl;
	output_matrix(error2);
}