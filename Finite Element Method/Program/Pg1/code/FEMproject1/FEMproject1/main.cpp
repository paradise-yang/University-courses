#include<iostream>
#include <iomanip>
#include<vector>
#include<cmath>

#include"femori.h"

using namespace std;

int main()
{
	vector<vector<double>> error(5, vector<double>(4, 0));
	error_out(error);
	cout << "基于一次有限元空间求解误差与收敛阶：" << endl;
	cout << "N   L无穷误差    误差收敛阶  L1误差       误差收敛阶" << endl;
	for (int i = 0; i < size(error); i++)
		cout << pow(2, i) * 10 << "  " << error[i][0] << "   " << error[i][1] << "    " << error[i][2] << "   " << error[i][3] << endl;

}