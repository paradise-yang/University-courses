#include<iostream>
#include<vector>
#include <cstdlib>
#include<cmath>
#include<iomanip>

using namespace std;

double f(double t, double x)
{
	return (t - exp(-t)) / (x + exp(x));
}

void adamsbashforth(vector<vector<double>>& num, double x, double y, int k)
{
	int N = pow(2, k);
	double h = 1.0 / N;
	vector<double> temp = { x, y };
	num.push_back(temp);
	//5阶Runge-Kutta计算前四项；
	for (int i = 1; i <= 4; i++)
	{
		double F1 = h * f(x, y);
		double F2 = h * f(x + h / 4, y + F1 / 4);
		double F3 = h * f(x + 3.0 / 8 * h, y + 3.0 / 32 * F1 + 9.0 / 32 * F2);
		double F4 = h * f(x + 12.0 / 13 * h, y + 1932.0 / 2197 * F1 - 7200.0 / 2197 * F2 + 7296.0 / 2197 * F3);
		double F5 = h * f(x + h, y + 439.0 / 216 * F1 - 8 * F2 + 3680.0 / 513 * F3 - 845.0 / 4104 * F4);
		double F6 = h * f(x + 1.0 / 2 * h, y - 8.0 / 27 * F1 + 2 * F2 - 3544.0 / 2565 * F3 + 1859.0 / 4104 * F4 - 11.0 / 40 * F5);
		double y5 = y + 16.0 / 135 * F1 + 6656.0 / 12825 * F3 + 28561.0 / 56430 * F4 - 9.0 / 50 * F5 + 2.0 / 55 * F6;
		x = h + x;
		y = y5;
		temp = { x, y };
		num.push_back(temp);
	}
	//Adams-Bashforth计算余下的项；
	for (int i = 5; i <= N; i++)
	{
		double tempy = num[i - 1][1] +
			h / 720 * 
			(1901 * f(num[i - 1][0], num[i - 1][1]) 
				- 2774 * f(num[i - 2][0], num[i - 2][1]) 
				+ 2616 * f(num[i - 3][0], num[i - 3][1]) 
				- 1274 * f(num[i - 4][0], num[i - 4][1]) 
				+ 251 * f(num[i - 5][0], num[i - 5][1]));
		double tempx = h + num[i - 1][0];
		temp = { tempx, tempy };
		num.push_back(temp);
	}
}

void order(vector<vector<double>>& ord)
{
	ord[0][1] = 0;
	for (int i = 1; i < ord.size(); i++)
		ord[i][1] = log(ord[i - 1][0] / ord[i][0]) / log(2);
}

int main()
{
	double  x = 0.0, y = 0.0;
	vector<vector<double>> ord(6, vector<double>(2, 0));
	for (int k = 3; k <= 8; k++)
	{
		vector<vector<double>> num;
		adamsbashforth(num, x, y, k);
		ord[k - 3][0] = fabs(-1.0- num[size(num) - 1][1]);
		//cout << pow(2, k) << "  " << ord[k - 3][0] << endl;
	}
	order(ord);
	cout << "当解x(1)=-1时相应的误差与阶为：" << endl << endl;
	cout << "  N       误差        阶" << endl;
	cout <<setw(3) << pow(2, 3) << "  " << setw(10) << ord[0][0] << "      - " << endl;
	for (int k = 4; k <= 8; k++)
		cout << setw(3) << pow(2, k) << "  " << setw(10) << ord[k - 3][0] << "  " << setw(10) << ord[k - 3][1] << endl;

	cout << endl << endl << endl;

	for (int k = 3; k <= 8; k++)
	{
		vector<vector<double>> num;
		adamsbashforth(num, x, y, k);
		ord[k - 3][0] = fabs(num[size(num) - 1][1]+ 0.1557824245884601849);
		//cout << pow(2, k) << "  " << ord[k - 3][0] << endl;
	}
	order(ord);
	cout << "当解x(1)=-0.1557824245时相应的误差与阶为：" << endl << endl;
	cout << "  N       误差        阶" << endl;
	cout << setw(3) << pow(2, 3) << "  " << setw(10) << ord[0][0] << "      - " << endl;
	for (int k = 4; k <= 8; k++)
		cout << setw(3) << pow(2, k) << "  " << setw(10) << ord[k - 3][0] << "  " << setw(10) << ord[k - 3][1] << endl;
}