#include<iostream>
#include <iomanip>
#include<vector>
#include<cmath>

using namespace std;

//注：若想查看具体数值解，请取消相关注释！

int main()
{
	//问题1.1。
	cout << "问题1.1";
	int J = 80;                                            //空间离散程度；
	double deltax = 1.0 / J;
	double lamuda = 0.5;                                   //固定\lamuda值；
	double deltat = lamuda / J;                            //求解空间离散程度；
	double pi = 3.1415926535897932384626433832795028841971;//圆周率\pi值；
	vector<vector<double>> u(2, vector<double>(J + 1, 0)); //记录求解向量；                     
	
	//FTCS格式
	cout << endl << "FTCS格式数值求解：" << endl;
	vector<double> error_t(10, 0);
	for (int i = 1; i <= 10; i++)
	{
		//初值；
		for (int j = 0; j <= J; j++)
			u[0][j] = sin(2 * pi * j * deltax);
		double T = i * 0.1;//待求某时间出数值解；
		int N = T * J / lamuda;//需要前进的步数；
		//数值格式求解；
		for (int i = 1; i <= N; i++)
		{
			u[i % 2][0] = u[(i + 1) % 2][0] + lamuda / 2 * (u[(i + 1) % 2][1] - u[(i + 1) % 2][J - 1]);
			for (int j = 1; j < J; j++)
				u[i % 2][j] = u[(i + 1) % 2][j] + lamuda / 2 * (u[(i + 1) % 2][j + 1] - u[(i + 1) % 2][j - 1]);
			u[i % 2][J] = u[(i + 1) % 2][J] + lamuda / 2 * (u[(i + 1) % 2][1] - u[(i + 1) % 2][J - 1]);
		}
		//求解模最大误差；
		double max = 0;
		for (int j = 0; j <= J; j++)
			if (fabs(sin(2 * pi * (j * deltax + T)) - u[N % 2][j]) > max)
				max = fabs(sin(2 * pi * (j * deltax + T)) - u[N % 2][j]);
		error_t[i - 1] = max;
		cout << "T=" << T << ", 最大误差为：" << error_t[i - 1] << endl;
		/*
		//输出数值解，绘图需要；
		cout << endl;
		if (T == 0.1 || T == 1.0 || T == 0.4 || T == 0.8)
			for (int j = 0; j <= J; j++)
				cout << "{" << j * deltax << "," << u[N % 2][j] << "},";
		cout << endl << endl;
		*/
	}
	//for (int j = 0; j < error_t.size(); j++)
	//	cout << "{" << 0.1*(j+1) << "," << error_t[j] << "},";

	//Lax-Friedrich格式
	cout << endl << "Lax-Friedrich格式数值求解：" << endl;
	vector<double> error_LF(10, 0);
	for (int i = 1; i <= 10; i++)
	{
		//初值；
		for (int j = 0; j <= J; j++)
			u[0][j] = sin(2 * pi * j * deltax);
		double T = i * 0.1;//待求某时间出数值解；
		int N = T * J / lamuda;//需要前进的步数；
		//数值格式求解；
		for (int i = 1; i <= N; i++)
		{
			u[i % 2][0] = (1 + lamuda) / 2 * u[(i + 1) % 2][1] + (1 - lamuda) / 2 * u[(i + 1) % 2][J - 1];
			for (int j = 1; j < J; j++)
				u[i % 2][j] = (1 + lamuda) / 2 * u[(i + 1) % 2][j + 1] + (1 - lamuda) / 2 * u[(i + 1) % 2][j - 1];
			u[i % 2][J] = (1 + lamuda) / 2 * u[(i + 1) % 2][1] + (1 - lamuda) / 2 * u[(i + 1) % 2][J - 1];
		}
		//求解模最大误差；
		double max = 0;
		for (int j = 0; j <= J; j++)
			if (fabs(sin(2 * pi * (j * deltax + T)) - u[N % 2][j]) > max)
				max = fabs(sin(2 * pi * (j * deltax + T)) - u[N % 2][j]);
		error_LF[i - 1] = max;
		cout << "T=" << T << ", 最大误差为：" << error_LF[i - 1] << endl;
	}

	//Lax-Wendroff格式
	cout << endl << "Lax-Wendroff格式数值求解：" << endl;
	vector<double> error_LW(10, 0);
	for (int i = 1; i <= 10; i++)
	{
		//初值；
		for (int j = 0; j <= J; j++)
			u[0][j] = sin(2 * pi * j * deltax);
		double T = i * 0.1;//待求某时间出数值解；
		int N = T * J / lamuda;//需要前进的步数；
		//数值格式求解；
		for (int i = 1; i <= N; i++)
		{
			u[i % 2][0] = (1 - lamuda * lamuda) * u[(i + 1) % 2][0] + (lamuda * (lamuda + 1) / 2) * u[(i + 1) % 2][1] + (lamuda * (lamuda - 1) / 2) * u[(i + 1) % 2][J - 1];
			for (int j = 1; j < J; j++)
				u[i % 2][j] = (1 - lamuda * lamuda) * u[(i + 1) % 2][j] + (lamuda * (lamuda + 1) / 2) * u[(i + 1) % 2][j + 1] + (lamuda * (lamuda - 1) / 2) * u[(i + 1) % 2][j - 1];
			u[i % 2][J] = (1 - lamuda * lamuda) * u[(i + 1) % 2][J] + (lamuda * (lamuda + 1) / 2) * u[(i + 1) % 2][1] + (lamuda * (lamuda - 1) / 2) * u[(i + 1) % 2][J - 1];
		}
		//求解模最大误差；
		double max = 0;
		for (int j = 0; j <= J; j++)
			if (fabs(sin(2 * pi * (j * deltax + T)) - u[N % 2][j]) > max)
				max = fabs(sin(2 * pi * (j * deltax + T)) - u[N % 2][j]);
		error_LW[i - 1] = max;
		cout << "T=" << T << ", 最大误差为：" << error_LW[i - 1] << endl;
	}


	//问题1.2；
	cout << endl << endl << "问题1.2";
	cout << endl << "Lax-Wendroff格式数值求解：" << endl;
	lamuda = 0.5;
	double T = 1.0;
	vector<double> error_LW_2(16, 0);
	for (int i = 0; i < 16; i++)
	{
		J = 10 * (i + 1);
		//J = pow(2, i) * 10;
		deltax = 1.0 / J;
		deltat = lamuda / J;
		vector<vector<double>> u2(2, vector<double>(J + 1, 0)); //记录求解向量；
		//初值；
		for (int j = 0; j <= J; j++)
			u2[0][j] = sin(2 * pi * j * deltax);
		int N = T * J / lamuda;//需要前进的步数；
		//数值格式求解；
		for (int i = 1; i <= N; i++)
		{
			u2[i % 2][0] = (1 - lamuda * lamuda) * u2[(i + 1) % 2][0] + (lamuda * (lamuda + 1) / 2) * u2[(i + 1) % 2][1] + (lamuda * (lamuda - 1) / 2) * u2[(i + 1) % 2][J - 1];
			for (int j = 1; j < J; j++)
				u2[i % 2][j] = (1 - lamuda * lamuda) * u2[(i + 1) % 2][j] + (lamuda * (lamuda + 1) / 2) * u2[(i + 1) % 2][j + 1] + (lamuda * (lamuda - 1) / 2) * u2[(i + 1) % 2][j - 1];
			u2[i % 2][J] = (1 - lamuda * lamuda) * u2[(i + 1) % 2][J] + (lamuda * (lamuda + 1) / 2) * u2[(i + 1) % 2][1] + (lamuda * (lamuda - 1) / 2) * u2[(i + 1) % 2][J - 1];
		}
		/*//输出数值解，绘图需要；
		cout << endl;
		if (J == 10 || J == 20 || J == 40 || J == 80 || J == 160)
			for (int j = 0; j <= J; j++)
				cout << "{" << j * deltax << "," << u2[N % 2][j] << "},";
		cout << endl;
		*/
		//求解模最大误差；
		double max = 0;
		for (int j = 0; j <= J; j++)
			if (fabs(sin(2 * pi * (j * deltax + T)) - u2[N % 2][j]) > max)
				max = fabs(sin(2 * pi * (j * deltax + T)) - u2[N % 2][j]);
		error_LW_2[i] = max;
		cout << "J=" << J << ", 最大误差为：" << error_LW_2[i] << endl;
	}
	//for (int j = 0; j < error_LW_2.size(); j++)
	//	cout << "{" << (1+j) * 10 << "," << error_LW_2[j] << "},";
}