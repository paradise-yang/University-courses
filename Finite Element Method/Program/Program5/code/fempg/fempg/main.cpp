#include"femori.h"

using namespace std;

int main()
{
	int n = 9;
	vector<vector<double>> error(n, vector<double>(6, 0));
	/*
	cout << "�����ݶȷ���⣺" << endl;
	cout << "N\t����ʱ��\t L�������\t���������\t L1���\t ���������" << endl;
	CG_error_out(error, n);
	matrix_output(error);
	*/

	cout << "����������⣺" << endl;
	cout << "N\t����ʱ��\t L�������\t���������\t L1���\t���������" << endl;
	vector<vector<double>> errors(n, vector<double>(6, 0));
	MG_error_out_two(errors, n);
	matrix_output(errors);

	/*
	cout << "��������two����⣺" << endl;
	cout << "N\t����ʱ��\t L�������\t���������\t L1���\t���������" << endl;
	vector<vector<double>> errorss(n, vector<double>(6, 0));
	MG_error_out_two(errorss, n);
	matrix_output(errorss);
	//*/
}