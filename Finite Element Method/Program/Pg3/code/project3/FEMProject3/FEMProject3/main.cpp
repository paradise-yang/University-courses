#include"FEMori.h"
#include"shinshkin.h"

using namespace std;

int main()
{
	int n = 7;
	vector<vector<double>> error(n, vector<double>(6, 0.0));
	error_out(error);
	/*
	cout << "����һ������Ԫ�ռ��������������ף�" << endl;
	cout << "N   \tL�������  \t���������  \tL1���  \t���������" << endl;
	for (int i = 0; i < size(error); i++)
		cout << pow(2, i) * 10 << "\t" << error[i][0] << "\t" << error[i][1] << "\t" << error[i][2] << "\t" << error[i][3] << endl;
	cout << endl;
	*/
	
	cout << "����shinshkin���ֵ�һ������Ԫ�ռ��������������ף�" << endl;
	cout << "N\t time\t L�������\t ���������\t L1���\t ���������" << endl;
	shinshkin_error_out(n, error);
	matrix_output(error);

}