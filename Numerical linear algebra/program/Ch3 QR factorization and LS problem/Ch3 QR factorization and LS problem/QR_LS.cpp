#include "QR_LS.h"

#define dim1 30   //Ch1.1��ά������֤����30��֮��ʼ��һ��ƫ�50�������ƫ��55���ٽ�ֵ����ά�� >55 ʱֻ��� -nan(ind)��
#define dim2 100  //Ch1.2.(1)��ά����
#define dim3 10   //Ch1.2.(2)��ά������֤����11���ٽ�ֵ����ά�� >11 ʱ�������Ѿ���ʼƫ��ϴ�

using namespace std;

typedef  vector<vector<double>> Matrix;

//��ע�������������㷨��
//һ������ʦ��˵�ġ�͵��������ʵ��Householder�����þ���˷�����QR�ֽ���LS���⣬
//��һ���ǲ�ֱ������Householder��ͨ��������������ϵ������QR�ֽ⣬
//��������ע��ԭ��ֻ�����˵�һ�ֽⷨ������뿴�ڶ���������������ͷ�ļ��и���ע�ͣ�����

int main()
{
	Matrix matrix, matrix2;

	//Ch3.1 QR�������Ch1�����������飬�������ý�����бȽϣ��Ƚ���ʵ�鱨��������������ﲻ���������Ch1�ļ���������
	cout << endl << "Ch3.1 QR�������Ch1�����������飬�������ý�����бȽ�." << endl;

	band_matrix1(dim1, matrix);//��ʼ������
	vector<double> b1(dim1);
	b1[0] = 7;                 //��ʼ������ֵ��
	b1[dim1-1] = 14;
	for (int i = 1; i < dim1-1; ++i)
		b1[i] = 15;
	LS(matrix, dim1, dim1, b1);
	cout << endl << "Ch1.1�����Է�����QR��ý��Ϊ��" << endl;
	Vector_cout(dim1,b1);
	matrix.clear();

	band_matrix2(dim2, matrix);//��ʼ������
	vector<double> b2(dim2);
	b2[0] = 11;                //��ʼ������ֵ��
	b2[dim2-1] = 11;
	for (int i = 1; i < dim2-1; ++i)
		b2[i] = 12;
	LS(matrix, dim2, dim2, b2);
	cout << endl << "Ch1.2.(1)�����Է�����QR��ý��Ϊ��" << endl;
	Vector_cout(dim2,b2);
	matrix.clear();

	Hilbert_matrix(dim3, matrix);
	vector<double> b3(dim3);
	for (int i = 0; i < dim3; i++)  //��ʼ������ֵ��
	{
		b3[i] = 0;
		for (int j = 0; j < dim3; j++)
			b3[i] += matrix[i][j];
	}
	LS(matrix, dim3, dim3, b3);
	cout << endl << "Ch1.2.(2)�����Է�����QR��ý��Ϊ��" << endl;
	Vector_cout(dim3,b3);
	cout << endl;
	matrix.clear();


	//Ch3.2 ����ζ���ʽ��ʹ�ò������ڵ�2-������С�����������Ӧ���ݣ�
	cout << endl << endl << "Ch3.2 ����ζ���ʽ��ʹ�ò�����2-������С�����������Ӧ����." << endl;
	vector<double> b4(7), y4(7);
	Ch3_2_make(matrix, b4);
	matrix2 = matrix; y4 = b4;
	LS(matrix, 7, 3, b4);
	cout << endl << "�������ʽΪ��y=" << b4[0] << " t^2 + " << b4[1] << " t + " << b4[2];
	cout << endl << "��ʱ������2-����Ϊ��||b - Ax|| = " << sqrt(norm_2(matrix2, 7, 3, y4, b4)) << endl;
	cout << endl;
	matrix.clear(); 
	matrix2.clear();


	//Ch3.3 �����������Ƶ�����ģ����С���˽⣻
	cout << endl << endl << "Ch3.3 �󷿲����Ƶ�����ģ����С���˽�." << endl;
	vector<double> b5(28), y5(28);
	Ch3_3_make(matrix, b5);
	matrix2 = matrix; y5 = b5;
	LS(matrix, 28, 12, b5);
	cout << endl << "������ģ�͵���С���˽�Ϊ��" << endl;
	Vector_cout(12, b5);
	cout << "��ʱ������2-����Ϊ��||b - Ax|| = " << sqrt(norm_2(matrix2, 28, 12, y5, b5)) << endl;
	cout << endl;
	matrix.clear();
	matrix2.clear();
}