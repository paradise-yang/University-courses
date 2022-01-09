#include<iostream>
#include <iomanip>
#include<vector>
#include<cmath>

using namespace std;

//注：若想查看具体数值解，请取消相关注释！

int main()
{
	//问题2.1。
	cout << "问题2.1" << endl;
	int J = 80;                                            //空间离散程度；
	double deltax = 1.0 / J;
	double T = 1.0;                                        //待求时间T处数值解；
	double lamuda = 0.5;                                   //\lamuda=deltat/deltax；
	int N = T * J / lamuda;                                //需要前进的步数；
	double deltat = lamuda / J;                            //求解空间离散程度；
	double pi = 3.1415926535897932384626433832795028841971;//圆周率\pi值；
	vector<vector<double>> u(3, vector<double>(J + 1, 0)); //记录求解向量；
	double max = 0;                                        //记录误差；
	
	//变化\lamuda。
	cout << "T=1.0,J=80,\lamuda=0.5:" << endl;
	lamuda = 0.5;
	J = 80;
	deltax = 1.0 / J;
	T = 1.0;
	N = T * J / lamuda;
	//初值；
	for (int j = 0; j <= J; j++)
		u[0][j] = sin(2 * pi * j * deltax);
	//第一步FTFS；
	for (int j = 0; j <= J; j++)
		u[1][j] = u[0][j] + lamuda / 2 * (u[0][(j + 1) % J] - u[0][j]);
	//CTCS格式；
	for (int i = 2; i <= N; i++)
	{
		u[i % 3][0] = u[(i - 2) % 3][0] + lamuda * (u[(i - 1) % 3][1] - u[(i - 1) % 3][J - 1]);
		for (int j = 1; j < J; j++)
			u[i % 3][j] = u[(i - 2) % 3][j] + lamuda * (u[(i - 1) % 3][j + 1] - u[(i - 1) % 3][j - 1]);
		u[i % 3][J] = u[(i - 2) % 3][J] + lamuda * (u[(i - 1) % 3][1] - u[(i - 1) % 3][J - 1]);
	}
	//求解模最大误差；
	max = 0;
	for (int j = 0; j <= J; j++)
		if (fabs(sin(2 * pi * (j * deltax + T)) - u[N % 3][j]) > max)
			max = fabs(sin(2 * pi * (j * deltax + T)) - u[N % 3][j]);
	cout << "模最大误差为：" << max << endl;
	//输出数值解；
	cout << "数值解结果：" << endl;
	for (int j = 0; j < size(u[0]); j++)
		cout << "{" << j * deltax << "," << u[N % 3][j] << "},";
	cout << endl << endl;

	cout << "T=1.0,J=80,\lamuda=1.5:" << endl;
	lamuda = 1.5;
	J = 80;
	deltax = 1.0 / J;
	T = 1.0;
	N = T * J / lamuda;
	double last_deltat = T - lamuda * deltax * N;
	N += 1;
	//初值；
	for (int j = 0; j <= J; j++)
		u[0][j] = sin(2 * pi * j * deltax);
	//第一步FTFS；
	for (int j = 0; j <= J; j++)
		u[1][j] = u[0][j] + lamuda / 2 * (u[0][(j + 1) % J] - u[0][j]);
	//CTCS格式；
	for (int i = 2; i < N; i++)
	{
		u[i % 3][0] = u[(i - 2) % 3][0] + lamuda * (u[(i - 1) % 3][1] - u[(i - 1) % 3][J - 1]);
		for (int j = 1; j < J; j++)
			u[i % 3][j] = u[(i - 2) % 3][j] + lamuda * (u[(i - 1) % 3][j + 1] - u[(i - 1) % 3][j - 1]);
		u[i % 3][J] = u[(i - 2) % 3][J] + lamuda * (u[(i - 1) % 3][1] - u[(i - 1) % 3][J - 1]);
	}
	//最后一步；
	u[N % 3][0] = u[(N - 2) % 3][0] + last_deltat / deltax * (u[(N - 1) % 3][1] - u[(N - 1) % 3][J - 1]);
	for (int j = 1; j < J; j++)
		u[N % 3][j] = u[(N - 2) % 3][j] + last_deltat / deltax * (u[(N - 1) % 3][j + 1] - u[(N - 1) % 3][j - 1]);
	u[N % 3][J] = u[(N - 2) % 3][J] + last_deltat / deltax * (u[(N - 1) % 3][1] - u[(N - 1) % 3][J - 1]);
	//求解模最大误差；
	max = 0;
	for (int j = 0; j <= J; j++)
		if (fabs(sin(2 * pi * (j * deltax + T)) - u[N % 3][j]) > max)
			max = fabs(sin(2 * pi * (j * deltax + T)) - u[N % 3][j]);
	cout << "模最大误差为：" << max << endl;
	//输出数值解；
	cout << "数值解结果：" << endl;
	for (int j = 0; j < size(u[0]); j++)
		cout << "{" << j * deltax << "," << u[N % 3][j] << "},";
	cout << endl << endl;



	//问题2.2。
	cout << endl << endl << endl << "问题2.2" << endl;
	cout << "T=1.0,J=80,\lamuda=0.5:" << endl;
	lamuda = 0.5;
	T = 1.0;
	for (int i = 0; i < 5; i++)
	{
		J = pow(2, i) * 10;
		deltax = 1.0 / J;
		deltat = lamuda / J;
		vector<vector<double>> u2(3, vector<double>(J + 1, 0)); //记录求解向量；
		//初值；
		for (int j = 0; j <= J; j++)
			u2[0][j] = sin(2 * pi * j * deltax);
		N = T * J / lamuda;
		//第一步FTFS；
		for (int j = 0; j <= J; j++)
			u2[1][j] = u2[0][j] + lamuda / 2 * (u2[0][(j + 1) % J] - u2[0][j]);
		//CTCS格式；
		for (int i = 2; i <= N; i++)
		{
			u2[i % 3][0] = u2[(i - 2) % 3][0] + lamuda * (u2[(i - 1) % 3][1] - u2[(i - 1) % 3][J - 1]);
			for (int j = 1; j < J; j++)
				u2[i % 3][j] = u2[(i - 2) % 3][j] + lamuda * (u2[(i - 1) % 3][j + 1] - u2[(i - 1) % 3][j - 1]);
			u2[i % 3][J] = u2[(i - 2) % 3][J] + lamuda * (u2[(i - 1) % 3][1] - u2[(i - 1) % 3][J - 1]);
		}
		//求解模最大误差；
		max = 0;
		for (int j = 0; j <= J; j++)
			if (fabs(sin(2 * pi * (j * deltax + T)) - u2[N % 3][j]) > max)
				max = fabs(sin(2 * pi * (j * deltax + T)) - u2[N % 3][j]);
		cout << "J=" << J << ",模最大误差为：" << max << endl;
		/*
		//输出数值解；
		cout << "数值解结果：" << endl;
		for (int j = 0; j < size(u2[0]); j++)
			cout << "{" << j * deltax << "," << u2[N % 3][j] << "},";
		cout << endl;
		*/
	}



	//问题2.3。
	cout << endl << endl << endl << endl << "问题2.3" << endl;
	cout << "\lamuda=0.5,J=80,T=0.2,FTBS结果:" << endl;
	J = 80;
	deltax = 1.0 / J;
	lamuda = 0.5;
	T = 0.2;
	N = T * J / lamuda;
	deltat = lamuda / J;
	vector<vector<double>> u3(2, vector<double>(J + 1, 0));
	//初值；
	for (int j = 0; j <= J; j++)
		u3[0][j] = sin(2 * pi * j * deltax);
	//FTBS格式；
	for (int i = 1; i <= N; i++)
	{
		u3[i % 2][0] = u3[(i + 1) % 2][0] + lamuda * (u3[(i + 1) % 2][0] - u3[(i + 1) % 2][J - 1]);
		for (int j = 1; j <= J; j++)
			u3[i % 2][j] = u3[(i + 1) % 2][j] + lamuda * (u3[(i + 1) % 2][j] - u3[(i + 1) % 2][j - 1]);
	}
	//求解模最大误差；
	max = 0;
	for (int j = 0; j <= J; j++)
		if (fabs(sin(2 * pi * (j * deltax + T)) - u3[N % 2][j]) > max)
			max = fabs(sin(2 * pi * (j * deltax + T)) - u3[N % 2][j]);
	cout << "模最大误差为：" << max << endl;
	//输出数值解；
	cout << "数值解结果：" << endl;
	for (int j = 0; j < size(u3[0]); j++)
		cout << "{" << j * deltax << "," << u3[N % 2][j] << "},";
	cout << endl << endl;

	cout << "\lamuda=0.5,J=80,T=0.5,FTBS结果:" << endl;
	J = 80;
	deltax = 1.0 / J;
	lamuda = 0.5;
	T = 0.5;
	N = T * J / lamuda;
	deltat = lamuda / J;
	//初值；
	for (int j = 0; j <= J; j++)
		u3[0][j] = sin(2 * pi * j * deltax);
	//FTBS格式；
	for (int i = 1; i <= N; i++)
	{
		u3[i % 2][0] = u3[(i + 1) % 2][0] + lamuda * (u3[(i + 1) % 2][0] - u3[(i + 1) % 2][J - 1]);
		for (int j = 1; j <= J; j++)
			u3[i % 2][j] = u3[(i + 1) % 2][j] + lamuda * (u3[(i + 1) % 2][j] - u3[(i + 1) % 2][j - 1]);
	}
	//求解模最大误差；
	max = 0;
	for (int j = 0; j <= J; j++)
		if (fabs(sin(2 * pi * (j * deltax + T)) - u3[N % 2][j]) > max)
			max = fabs(sin(2 * pi * (j * deltax + T)) - u3[N % 2][j]);
	cout << "模最大误差为：" << max << endl;
	//输出数值解；
	cout << "数值解结果：" << endl;
	for (int j = 0; j < size(u3[0]); j++)
		cout << "{" << j * deltax << "," << u3[N % 2][j] << "},";
	cout << endl << endl;
}