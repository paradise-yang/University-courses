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
	cout << "����һ������Ԫ�ռ��������������ף�" << endl;
	cout << "N   L�������    ���������  L1���       ���������" << endl;
	for (int i = 0; i < size(error); i++)
		cout << pow(2, i) * 10 << "  " << error[i][0] << "   " << error[i][1] << "    " << error[i][2] << "   " << error[i][3] << endl;

}