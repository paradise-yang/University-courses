#include "Symmetric_eigenvalue.h"

using namespace std;

typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

int main()
{
	//Ch7.1��Jacobi���ط���50-100�׶Գ��������ȫ������ֵ��
	cout << "Ch7.1 ��Jacobi���ط���50-100�׶Գ��������ȫ������ֵ��" << endl << endl;
	cout << "���ʣ������ĳ��������������ֵ������50-100�����о������ ǰ������1����������0��";
	int temp = 0;
	double u = 1e-7;//���ȣ�Ҳ����ֵdeta��
	cin >> temp;
	if (temp == 1)
	{
		int dim = 0;
		cout << endl << "��������鿴�ľ��������";
		cin >> dim;
		Matrix matrix(dim, Vector(dim, 0));
		for (int i = 0; i < dim; i++) matrix[i][i] = 4;
		for (int i = 1; i < dim; i++) matrix[i][i - 1] = 1, matrix[i - 1][i] = 1;
		Pass_Jacobi_method(matrix, u);
		Matrix_cout(matrix, dim, dim);
	}
	else
	{
		for (int dim = 50; dim <= 100; dim++)
		{
			Matrix matrix(dim, Vector(dim, 0));
			for (int i = 0; i < dim; i++) matrix[i][i] = 4;
			for (int i = 1; i < dim; i++) matrix[i][i - 1] = 1, matrix[i - 1][i] = 1;
			Pass_Jacobi_method(matrix, u);
		}
	}



	//Ch7.2�ö��ַ���ʵ�Գ�������ָ������ֵ��
	cout << endl << endl << endl;
	cout << "Ch7.2 �ö��ַ���ʵ�Գ�������ָ������ֵ��" << endl;
	int dim = 100;
	Matrix matrix(dim, Vector(dim, 0));
	for (int i = 0; i < dim; i++) matrix[i][i] = 2;
	for (int i = 1; i < dim; i++) matrix[i][i - 1] = -1, matrix[i - 1][i] = -1;
	cout << "�����һ������ֵ����С����ֵ��Ϊ " << bisection_method(matrix, 1) << endl;
	cout << "�����100������ֵ���������ֵ��Ϊ " << bisection_method(matrix, dim) << endl;
}