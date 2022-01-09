#include<iostream>
#include <iomanip>
#include<vector>
#include<cmath>

using namespace std;

int main()
{
	int J = 50;
	double deltax = 0.02;
	int N = 30;
	double deltat = 0.01;
	double pi = 3.141592653589793238;
	vector<vector<double>> u(2, vector<double>(J + 1, 0));

	//问题1
	//边界初值；
	for (int j = 0; j <= J; j++)
		u[0][j] = sin(2 * pi * j * deltax);
	//差分；
	for (int i = 1; i <= N; i++)
	{
		for (int j = 0; j < J; j++)
			u[i % 2][j % J] = u[(i + 1) % 2][j % J] + deltat / deltax * (u[(i + 1) % 2][(j + 1) % J] - u[(i + 1) % 2][j % J]);
		u[i % 2][J] = u[(i + 1) % 2][J] + deltat / deltax * (u[(i + 1) % 2][1] - u[(i + 1) % 2][J]);
	}
	//输出；
	cout << "deltat=0.01时(x坐标，差分解，精确解）" << endl;
	for (int j = 0; j <= J; j++)
	{
		cout << "(" << setprecision(3) << j * deltax << ",   " << setprecision(6) << u[0][j] << "，   " << setprecision(6) << sin(2 * pi * (j * deltax + 0.3)) << ")     ";
		if (j % 3 == 2) cout << endl;
	}
	//绘图所需输出
	//for (int j = 0; j <= J; j++)
	//	cout << "{" << j * deltax << "," << u[0][j] << "},";

	
	cout << endl << endl;


	//问题2；
	deltat = 0.03;
	N = 10;
	//边界初值；
	for (int j = 0; j <= J; j++)
		u[0][j] = sin(2 * pi * j * deltax);
	//差分；
	for (int i = 1; i <= N; i++)
	{
		for (int j = 0; j < J; j++)
			u[i % 2][j % J] = u[(i + 1) % 2][j % J] + deltat / deltax * (u[(i + 1) % 2][(j + 1) % J] - u[(i + 1) % 2][j % J]);
		u[i % 2][J] = u[(i + 1) % 2][J] + deltat / deltax * (u[(i + 1) % 2][1] - u[(i + 1) % 2][J]);
	}
	//输出；
	cout << "deltat=0.03时(x坐标，差分解，精确解）" << endl;
	for (int j = 0; j <= J; j++)
	{
		cout << "(" << setprecision(3) << j * deltax << ",   " << setprecision(6) << u[0][j] << "，   " << setprecision(6) << sin(2 * pi * (j * deltax + 0.3)) << ")     ";
		if (j % 3 == 2) cout << endl;
	}
	//绘图所需输出
	//for (int j = 0; j <= J; j++)
	//	cout << "{" << j * deltax << "," << u[0][j] << "},";

}