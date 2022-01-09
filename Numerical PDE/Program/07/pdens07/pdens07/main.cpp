#include<iostream>
#include<iomanip>
#include<vector>
#include<cmath>
#include<fstream>
#include<string>

using namespace std;

#define pi 3.1415926535897932384626433832795028841971 //圆周率\pi值；
#define miu 0.1 //网格比；

//v(x,0)=f(x)=x(1-x)
double f(double x){	return x * (1 - x);}
//v_x(0,t)=a(t)=10sin t
double a(double t){	return 10.0*sin(t);}
//v(1,t)=b(t)=4sin 6t
double b(double t){	return 4.0*sin(6.0*t);}
//F(x,t)=sin 2pi x * sin 4pi t
double F(double x, double t){ return sin(2 * pi*x)*sin(4 * pi*t);}

void FTCS_first(double v, int J, double deltat, double T, vector<vector<double>> &u)
{
	double deltax = 1.0 / J;
	double lambda = v * deltat / deltax / deltax;
	int N = T / deltat;
	//初值；
	for (int j = 0; j <= J; j++)
		u[0][j] = f(j * deltax);
	//FTCS格式；
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j < J; j++)
			u[i % 2][j] = u[(i + 1) % 2][j] + lambda * (u[(i + 1) % 2][j + 1] - 2 * u[(i + 1) % 2][j] + u[(i + 1) % 2][j - 1]) + deltat * F(j*deltax, i*deltat);
		//左边界处理；
		u[i % 2][0] = u[i % 2][1] - deltax * a(i*deltat);
		//右边界处理；
		u[i % 2][J] = b(i*deltat);
	}
}

void FTCS_second(double v, int J, double deltat, double T, vector<vector<double>> &u)
{
	double deltax = 1.0 / J;
	double lambda = v * deltat / deltax / deltax;
	int N = T / deltat;
	//初值；
	for (int j = 0; j <= J; j++)
		u[0][j] = f(j * deltax);
	//FTCS格式；
	for (int i = 1; i <= N; i++)
	{
		for (int j = 1; j < J; j++)
			u[i % 2][j] = u[(i + 1) % 2][j] + lambda * (u[(i + 1) % 2][j + 1] - 2 * u[(i + 1) % 2][j] + u[(i + 1) % 2][j - 1]) + deltat * F(j*deltax, i*deltat);
		//左边界处理；
		u[i % 2][0] = u[(i + 1) % 2][0] + 2 * lambda*(u[(i + 1) % 2][1] - u[(i + 1) % 2][0]) - 2 * lambda*deltax* a((i - 1)*deltat);// -deltax * a(i*deltat);
		//右边界处理；
		u[i % 2][J] = b(i*deltat);
	}
}

//v222(x,0)=sin 4pi x
double v222(double x) { return sin(4 * pi*x);}
void FTCS_third(double v, double deltax, double deltat, double T, vector<vector<double>> &u)
{
	int J = 1.0 / deltax;
	double lambda = v * deltat / deltax / deltax;
	int N = T / deltat;
	//初值；
	for (int j = 0; j <= J; j++)
		u[0][j] = v222(j * deltax);
	//FTCS格式；
	for (int i = 1; i <= N; i++)
	{
		//左边界处理；
		u[i % 2][0] = 0;
		for (int j = 1; j < J; j++)
			u[i % 2][j] = u[(i + 1) % 2][j] + lambda * (u[(i + 1) % 2][j + 1] - 2 * u[(i + 1) % 2][j] + u[(i + 1) % 2][j - 1]);
		//右边界处理；
		u[i % 2][J] = 0;
	}
}

//输出
void output_file(int i, int sign, int J, double deltat, double T, vector<vector<double>> &u)
{
	double deltax = 1.0 / J;
	int N = T / deltat;
	ofstream outfile;
	string temp = "E:\\study_materials\\PDEns\\07\\";
	string path;
	if (sign == 1)
		path = temp + "first" + to_string(i) + ".txt";
	else
		path = temp + "second" + to_string(i) + ".txt";
	outfile.open(path, ios::out);
	for (int j = 0; j < size(u[0]); j++)
		outfile << j * deltax << " " << u[N % 2][j] << endl;
	outfile.close();
}

int main()
{
	//HW1.5.9;
	double v = 0.1;
	int J = 10;
	vector<double> t = { 0.1,0.9,2.0 };
	double deltat = 0.05;
	vector<vector<double>> u_first(2, vector<double>(J + 1, 0)), u_second(2, vector<double>(J + 1, 0));
	
	for (int i = 0; i < 3; i++)
	{
		FTCS_first(v, J, deltat, t[i], u_first);
		output_file(i + 1, 1, J, deltat, t[i], u_first);
		FTCS_second(v, J, deltat, t[i], u_second);
		output_file(i + 1, 2, J, deltat, t[i], u_second);
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//HW1.5.10;
	J = 20;
	deltat = 0.01;
	vector<vector<double>> temp(2, vector<double>(J + 1, 0));
	u_first.clear(), u_second.clear();
	u_first = temp, u_second = temp;
	for (int i = 0; i < 3; i++)
	{
		FTCS_first(v, J, deltat, t[i], u_first);
		output_file(i + 4, 1, J, deltat, t[i], u_first);
		FTCS_second(v, J, deltat, t[i], u_second);
		output_file(i + 4, 2, J, deltat, t[i], u_second);
	}
	
	J = 40;
	deltat = 0.002;
	vector<vector<double>> tempp(2, vector<double>(J + 1, 0));
	u_first.clear(), u_second.clear();
	u_first = tempp, u_second = tempp;
	for (int i = 0; i < 3; i++)
	{
		FTCS_first(v, J, deltat, t[i], u_first);
		output_file(i + 7, 1, J, deltat, t[i], u_first);
		FTCS_second(v, J, deltat, t[i], u_second);
		output_file(i + 7, 2, J, deltat, t[i], u_second);
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//HW2.2.2;
	//(i)
	t[0] = 0.05, t[1] = 0.1;
	double deltax = 0.1;
	J = 1.0 / deltax;
	deltat = 0.05;
	vector<vector<double>> temppp(2, vector<double>(J + 1, 0));
	u_first.clear(), u_first = temppp;
	for (int i = 0; i < 2; i++)
		FTCS_third(v, deltax, deltat, t[i], u_first),
		output_file(i + 10, 1, J, deltat, t[i], u_first);

	//(ii)
	deltax = 0.05, deltat = 0.0125, J = 1.0 / deltax;
	vector<vector<double>> tempppp(2, vector<double>(J + 1, 0));
	u_first.clear(), u_first = tempppp;
	for (int i = 0; i < 2; i++)
		FTCS_third(v, deltax, deltat, t[i], u_first),
		output_file(i + 12, 1, J, deltat, t[i], u_first);

	//(iii)
	deltax = 0.01, deltat = 0.0005, J = 1.0 / deltax;
	vector<vector<double>> temppppp(2, vector<double>(J + 1, 0));
	u_first.clear(), u_first = temppppp;
	for (int i = 0; i < 2; i++)
		FTCS_third(v, deltax, deltat, t[i], u_first),
		output_file(i + 14, 1, J, deltat, t[i], u_first);

	cout << "输出参见draft.nb文件" << endl;
}