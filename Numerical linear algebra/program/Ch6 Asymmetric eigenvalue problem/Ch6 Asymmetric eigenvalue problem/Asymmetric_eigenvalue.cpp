#include "Asymmetric_eigenvalue.h"

using namespace std;

typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

int main()
{
	//Ch6.1Ӧ���ݷ�������ʽģ������
	cout << "Ch6.1 Ӧ���ݷ�������ʽģ������" << endl;
	double error = 1e-6;
	Vector poly1 = { 1,1,-5,3 };
	cout << "The maxmium modulus root of the first polynomial is " << max_modulus_root_polynomial(poly1, error) << endl;
	Vector poly2 = { 1,0,-3,-1 };
	cout << "The maxmium modulus root of the second polynomial is " << max_modulus_root_polynomial(poly2, error) << endl;
	Vector poly3 = { 1,101,208.01,10891.01,9802.08,79108.9,-99902,790,-1000 };
	cout << "The maxmium modulus root of the third polynomial is " << max_modulus_root_polynomial(poly3, error) << endl;



	//Ch6.2Ӧ��QR�㷨���ʵ�������е�����ֵ��

	//Ch6.2.2��ⷽ�����еĸ���
	cout << endl << endl << endl << "Ch6.2.2 ��ⷽ��x^41 + x^3 +1 = 0ȫ������" << endl;
	cout << "�÷���ȫ����Ϊ��" << endl;
	double u = 1.0e-10;
	Matrix matrix0(41, Vector(41, 0));
	for (int i = 1; i < 41; i++) 
		matrix0[i][i - 1] = 1;
	matrix0[0][40] = -1, matrix0[3][40] = -1;
	Real_Schur_decomposition(matrix0, u);
	Eigenvalue(matrix0);
	//Matrix_cout(matrix0, 41, 41);//����任���41��S���뷽��
	
	//Ch6.2.3����������ֵ�����۲�����ֵ�仯��
	cout << endl << endl << endl << "Ch6.2.3����������ֵ�����۲�����ֵ�仯��" << endl;
	cout << endl << "x = 0.9ʱ��" << endl;
	Matrix matrix = { {9.1,3.0,2.6,4.0},{4.2,5.3,4.7,1.6},{3.2,1.7,9.4,0.9},{6.1,4.9,3.5,6.2} };
	Real_Schur_decomposition(matrix, u);
	cout << "�����Schur�ֽ�Ϊ��" << endl;
	Matrix_cout(matrix, 4, 4);
	cout << "������������ֵΪ��";
	Eigenvalue(matrix);

	cout << endl << "x = 1.0ʱ��" << endl;
	matrix = { {9.1,3.0,2.6,4.0},{4.2,5.3,4.7,1.6},{3.2,1.7,9.4,1.0},{6.1,4.9,3.5,6.2} };
	Real_Schur_decomposition(matrix, u);
	cout << "�����Schur�ֽ�Ϊ��" << endl;
	Matrix_cout(matrix, 4, 4);
	cout << "������������ֵΪ��";
	Eigenvalue(matrix);

	cout << endl << "x = 1.1ʱ��" << endl;
	matrix = { {9.1,3.0,2.6,4.0},{4.2,5.3,4.7,1.6},{3.2,1.7,9.4,1.1},{6.1,4.9,3.5,6.2} };
	Real_Schur_decomposition(matrix, u);
	cout << "�����Schur�ֽ�Ϊ��" << endl;
	Matrix_cout(matrix, 4, 4);
	cout << "������������ֵΪ��";
	Eigenvalue(matrix);

	cout << endl << endl << "С̷test is sb:" << endl;
	Matrix temptest = { {6,5,-5},{2,6,-2},{2,5,-1} };
	double value = 0;
	Power_method(temptest.size(), temptest, 1e-8, value);
	cout << value << endl;
}