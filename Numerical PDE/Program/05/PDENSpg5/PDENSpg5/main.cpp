#include<iostream>
#include<iomanip>
#include<vector>
#include<cmath>
#include<fstream>
#include<string>

using namespace std;

#define pi 3.1415926535897932384626433832795028841971 //圆周率\pi值；

//初值函数；
double u_0(double x)
{
	if (x >= 0.4 && x <= 0.6)
		return 1.0;
	else
		return 0.0;
}
//精确解函数；
double u_exact(double x, double t)
{
	if ((x - t >= 0.4) && (x - t <= 0.6))
		return 1.0;
	else
		return 0.0;
}
//FTBS格式；
void FTBS(double deltax, double lamuda, double T, vector<vector<double>> &u)
{
	int J = 7.0 / deltax;           //空间所需网格数；
	double deltat = lamuda * deltax;//求解空间离散程度；
	int N = T / deltat;             //需要前进的步数；
	//初值；
	for (int j = 0; j <= J; j++)
		u[0][j] = u_0(j*deltax);
	//FTBS格式；
	for (int i = 1; i <= N; i++)
	{
		u[i % 2][0] = u[(i + 1) % 2][0] - lamuda * (u[(i + 1) % 2][0] - 0);//u[(i + 1) % 2][J - 1]);
		for (int j = 1; j <= J; j++)
			u[i % 2][j] = u[(i + 1) % 2][j] - lamuda * (u[(i + 1) % 2][j] - u[(i + 1) % 2][j - 1]);
	}
}
//求解模最大误差；
double error_max(double deltax, double lamuda, double T, vector<vector<double>> u)
{
	int J = 6.0 / deltax;           //空间所需网格数；
	double deltat = lamuda * deltax;//求解空间离散程度；
	int N = T / deltat;             //需要前进的步数；
	double max = 0;                 //记录误差；
	for (int j = 0; j <= J; j++)
		if (fabs(u_exact(j*deltax, T) - u[N % 2][j]) > max)
			max = fabs(sin(2 * pi * (j * deltax + T)) - u[N % 2][j]);
	return max;
}
//输出数值解进入指定文本文件；
void output(int i, double deltax, double lamuda, double T, vector<vector<double>> u)
{
	double deltat = lamuda * deltax;//求解空间离散程度；
	int N = T / deltat;             //需要前进的步数；
	ofstream outfile;//输出三维数值误差；
	string temp = "E:\\study_materials\\PDEns\\05\\ns";
	string path = temp + to_string(i) + ".txt";
	outfile.open(path, ios::out);
	for (int j = 0; j < size(u[0]); j++)
		outfile << j * deltax << " " << u[N % 2][j] << endl;
}

int main()
{
	double deltax = 0.05;                                  //空间离散程度；
	int J = 7.0 / deltax;
	double T = 1.0;                                        //待求时间T处数值解；
	double lamuda = 0.2;                                   //r=\lamuda=deltat/deltax；
	vector<vector<double>> u(3, vector<double>(J + 1, 0)); //记录求解向量；
	double max = 0;                                        //记录误差；
	
	cout << endl << "deltax=0.05,lamuda=0.2,T=1.0：" << endl;
	FTBS(deltax, lamuda, T, u);
	max = error_max(deltax, lamuda, T, u);
	cout << "模最大误差为：" << max << endl;
	//输出数值解；
	output(0, deltax, lamuda, T, u);
	cout << "数值解结果参见文档" << endl << endl;
	
	cout << "deltax=0.05,lamuda=0.2,T=2.0：" << endl;
	T = 2.0;
	FTBS(deltax, lamuda, T, u);
	max = error_max(deltax, lamuda, T, u);
	cout << "模最大误差为：" << max << endl;
	//输出数值解；
	output(1, deltax, lamuda, T, u);
	cout << "数值解结果参见文档" << endl << endl;

	cout << "deltax=0.05,lamuda=0.2,T=5.0：" << endl;
	T = 5.0;
	FTBS(deltax, lamuda, T, u);
	max = error_max(deltax, lamuda, T, u);
	cout << "模最大误差为：" << max << endl;
	//输出数值解；
	output(2, deltax, lamuda, T, u);
	cout << "数值解结果参见文档" << endl << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	cout << endl << "deltax=0.05,lamuda=0.8,T=1.0：" << endl;
	lamuda = 0.8;
	T = 1.0;
	FTBS(deltax, lamuda, T, u);
	max = error_max(deltax, lamuda, T, u);
	cout << "模最大误差为：" << max << endl;
	//输出数值解；
	output(3, deltax, lamuda, T, u);
	cout << "数值解结果参见文档" << endl << endl;

	cout << "deltax=0.05,lamuda=0.8,T=2.0：" << endl;
	T = 2.0;
	FTBS(deltax, lamuda, T, u);
	max = error_max(deltax, lamuda, T, u);
	cout << "模最大误差为：" << max << endl;
	//输出数值解；
	output(4, deltax, lamuda, T, u);
	cout << "数值解结果参见文档" << endl << endl;

	cout << "deltax=0.05,lamuda=0.8,T=5.0：" << endl;
	T = 5.0;
	FTBS(deltax, lamuda, T, u);
	max = error_max(deltax, lamuda, T, u);
	cout << "模最大误差为：" << max << endl;
	//输出数值解；
	output(5, deltax, lamuda, T, u);
	cout << "数值解结果参见文档" << endl << endl;
}