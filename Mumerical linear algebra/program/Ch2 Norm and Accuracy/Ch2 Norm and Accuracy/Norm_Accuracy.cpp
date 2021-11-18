#include "Norm_Accuracy.h"

using namespace std;

typedef  vector<vector<double>> Matrix;

int main()
{
	Matrix matrix;

	//Chapter2.1 ����5-20��Hilbert������������
	cout << endl;
	cout << "Chapter2.1 ����5-20��Hilbert����������������" << endl;
	cout << endl;
	for (int n = 5; n <= 20; n++)
	{
		Hilbert_matrix(n, matrix); //����n��Hilbert����
		cout << n << "���������������" << Conditional_number_infinity(matrix) << endl;
		matrix.clear();
	}
	cout << endl;

	//Chapter2.2 ���ѡȡx,����b=Ax,Ȼ��������ԪGuass��ȥ����ⷽ���飬����5-30�׽�ľ��Ȳ�����ʵ���������Ƚϣ�
	//���ú���ѡȡ��xα�����������ѡȡx��Ĭ��Ϊ2^(1/2)*i����θ���ע�ͺ󣬿���ȡ2^(1/2)rand()��
	cout << endl;
	cout << "Chapter2.2 ���ѡȡx,����b=Ax,Ȼ��������ԪGuass��ȥ����ⷽ���飬����5-30�׽�ľ��Ȳ�����ʵ���������Ƚ�" << endl;
	cout << endl;
	for (int n = 5; n <= 30; n++)
	{
		vector<double> x(n), b(n), value(n),y(n);
		Vector_creat_random(n, x);                 //���������x��
		creat_matrix(n, matrix);                   //����ϵ������A��
		cout << n << "�������ѡȡ��xΪ��" << endl;
		Vector_cout(x);
		cout << endl;
		Matrix_times_vector(matrix, x, b);         //����b��
		Guass_Line_solution(matrix, b, value);     //���㷽�̵Ľ⣻
		Matrix_times_vector(matrix, value, y);     //����A���Լ����x��
		Vector_subtraction(y, b);                  //����b-Ax��������y���棻
		double accuracy = Conditional_number_infinity(matrix) * vector_infinity_norm(y) / vector_infinity_norm(b);
		cout << n << "�׽�ľ���Ϊ��" << accuracy << endl;
		Vector_subtraction(value,x);               //���㾫ȷ��x������Ĳ�ֵ��
		double real_accuracy = vector_infinity_norm(value) / vector_infinity_norm(x);
		cout << n << "��ʵ�����Ϊ��" << real_accuracy << endl << endl << endl;
		matrix.clear();
	}
}