#include<iostream>
#include<vector>
#include <cstdlib>
#include<cmath>

using namespace std;

double f(double x, double y)
{
	return exp(x * y) + cos(y - x);
}

void adapt(vector<vector<double>>& num, double x, double y, double delta, double h)
{
	//while(1)
	for(int i=0;i<=1000;i++)
	{
		double F1 = h*f(x, y);
		double F2 = h*f(x + h / 4, y + F1 / 4);
		double F3 = h * f(x + 3.0 / 8 * h, y + 3.0 / 32 * F1 + 9.0 / 32 * F2);
		double F4 = h * f(x + 12.0 / 13 * h, y + 1932.0 / 2197 * F1 - 7200.0 / 2197 * F2 + 7296.0 / 2197 * F3);
		double F5 = h * f(x + h, y + 439.0 / 216 * F1 - 8 * F2 + 3680.0 / 513 * F3 - 845.0 / 4104 * F4);
		double F6 = h * f(x + 1.0 / 2 * h, y - 8.0 / 27 * F1 + 2 * F2 - 3544.0 / 2565 * F3 + 1859.0 / 4104 * F4 - 11.0 / 40 * F5);
		double y5 = y + 16.0 / 135 * F1 + 6656.0 / 12825 * F3 + 28561.0 / 56430 * F4 - 9.0 / 50 * F5 + 2.0 / 55 * F6;
		double y4 = y + 25.0 / 216 * F1 + 1408.0 / 2565 * F3 + 2197.0 / 4104 * F4 - F5 / 5;
		if (F1 > 600|| F2 > 600|| F3 > 600)
			break;
		x = h + x;
		y = y5;
		vector<double> temp = { x, y };
		num.push_back(temp);
		double e = fabs(y5 - y4);
		h = 0.9 * h * pow(delta / e, 0.2);
	}
}

int main()
{
	double delta = 1e-7;
	double h = 0.01;
	double x = 1.0, y = 3.0;
	vector<vector<double>> num;
	adapt(num, x, y, delta, h);
	for (int i = 0; i < num.size(); i++)
	{
		cout << num[i][0] << " " << num[i][1] << endl;
	}
}