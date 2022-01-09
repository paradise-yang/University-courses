#include"uniform.h"
#include"shinshkin.h"

using namespace std;

int main()
{
	int n = 8;
	vector<vector<double>> error(n, vector<double>(8, 0.0));
	/* 
	cout << "\epsilon=0.1����������" << endl << endl;

	cout << "����һ������Ԫ�ռ��������������ף�" << endl;
	cout << "N\t time\t L�������\t ���������\t L1���\t ���������\t L2���\t���������" << endl;
	CG_error_out(error, n, epsilon1, 1);
	matrix_output(error);
	cout << "���ڶ�������Ԫ�ռ��������������ף�" << endl;
	cout << "N\t time\t L�������\t ���������\t L1���\t ���������\t L2���\t���������" << endl;
	CG_error_out_2(error, n, epsilon1, 1);
	matrix_output(error);
	cout << endl;
	cout << "*******************************************************************************************************************" << endl;
	cout << endl << endl;
	//*/
	/* 
	cout << "\epsilon=1e-7����������" << endl << endl;

	cout << "����һ������Ԫ�ռ��������������ף�" << endl;
	cout << "N\t time\t L�������\t ���������\t L1���\t ���������\t L2���\t���������" << endl;
	CG_error_out(error, n, epsilon7, 7);
	matrix_output(error);
	cout << "���ڶ�������Ԫ�ռ��������������ף�" << endl;
	cout << "N\t time\t L�������\t ���������\t L1���\t ���������\t L2���\t���������" << endl;
	CG_error_out_2(error, n, epsilon7, 7);
	matrix_output(error);
	cout << endl;
	cout << "*******************************************************************************************************************" << endl;
	cout << endl << endl;
	//*/
	///* 
	cout << "\epsilon=1e-7��shinshkin����" << endl << endl;

	cout << "����shinshkin���ֵ�һ������Ԫ�ռ��������������ף�" << endl;
	cout << "N\t time\t L�������\t ���������\t L1���\t ���������" << endl;
	shinshkin_error_out(n, error, 7);
	matrix_output(error);
	cout << "����shinshkin���ֵĶ�������Ԫ�ռ��������������ף�" << endl;
	cout << "N\t time\t L�������\t ���������\t L1���\t ���������" << endl;
	shinshkin_error_out_2(n, error, 7);
	matrix_output(error);
	cout << endl;
	//*/
}