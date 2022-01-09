#include "Conjugate_gradient_method.h"

using namespace std;

typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

int main() 
{
	int dim = 20;
	int times = 0;
	double w = 1.725;            //SOR�������ɳ����ӣ�
	double error = 1e-7;         //SOR�����ľ��ȣ�
	double epsilon = 1e-10;      //ʵ�ù����ݶȵ�������epsiloon��
	double iteration_error = 0;  //�����ݶȵ�����������ο���
	clock_t start_time, end_time;
	Matrix matrix((dim - 1)* (dim - 1), Vector((dim - 1)* (dim - 1), 0));
	Vector b((dim - 1) * (dim - 1), 0), x((dim - 1) * (dim - 1), 0);
	

	//Ch5.1 �����ݶȷ����PDE������SOR�����Ƚϣ�Ѱ������ɳ����ӣ�
	cout << "Ch5.1 �����ݶȷ����PDE������SOR�����Ƚϣ�Ѱ������ɳ����ӡ�" << endl;
	Ch5_1(dim, matrix, b);
	Conjugate_Gradient((dim - 1) * (dim - 1), matrix, b, epsilon, x, times, iteration_error);
	cout << endl << "�����ݶȷ���������Ϊ = " << times << endl;
	cout << "�ο����Ϊ = " << iteration_error << endl;
	cout << "������Ϊ��" << endl;
	Vector_cout(x);
	
	cout << endl;
	Vector x0((dim - 1) * (dim - 1), 1.0);
	start_time = clock();                //��ȡ��һ��ʱ�䣻
	SOR_iteration((dim - 1) * (dim - 1), w, matrix, x0, b, error, times);
	//Vector_cout(x0);
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "SOR�������ɳ�����Ϊw = " << w << endl;
	cout << "��������Ϊ = " << times << endl;
	cout << "����ʱ��Ϊ =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	


	//Ch5.2 ��Hilbert������Թ����ݶȷ���
	cout << endl << endl <<endl << "Ch5.2 ��Hilbert������Թ����ݶȷ���" << endl;
	
	dim = 100;     //Hilbert����ά����
	Matrix hilbert(dim,Vector(dim,0));
	Vector value(dim, 0), y(dim, 0);
	Hilbert_matrix(dim, hilbert);
	for (int i = 0; i < dim; i++)                //��ʼ��ֵ��
	{
		for (int j = 0; j < dim; j++)
			value[i] += hilbert[i][j];
		value[i] = value[i] / 3.0;
	}
	start_time = clock();                //��ȡ��һ��ʱ�䣻
	Conjugate_Gradient(dim, hilbert, value, epsilon, y, times, iteration_error);
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << endl << dim << "��Hilbert����" << "����������" << times << endl;
	cout << "����ʱ��Ϊ =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	cout << "�ο����Ϊ = " << iteration_error << endl;
	cout << "������Ϊ��" << endl;
	Vector_cout(y);
	


	//Ch5.3�ֱ���Jacobi��Gauss-Seidel�������ݶȷ���ⷽ���飻
	cout << endl << endl << endl << "Ch5.3�ֱ���Jacobi��Gauss-Seidel�������ݶȷ���ⷽ���顣" << endl;
	dim = 5;
	Matrix A = { {10.0,1.0,2.0,3.0,4.0},{1.0,9.0,-1.0,2.0,-3.0},{2.0,-1.0,7.0,3.0,-5.0},{3.0,2.0,3.0,12.0,-1.0},{4.0,-3.0,-5.0,-1.0,15.0} };
	Vector z = { 12.0,-27.0,14.0,-17.0,12.0 };
	Vector z0(dim, 1.0);
	
	cout << endl << "Jacobi������" << endl;
	Jacobi_iteration(dim, A, z0, z, error, times);
	cout << "��������Ϊ = " << times << endl;
	cout << "������Ϊ��";
	Vector_cout(z0);
	
	A = { {10.0,1.0,2.0,3.0,4.0},{1.0,9.0,-1.0,2.0,-3.0},{2.0,-1.0,7.0,3.0,-5.0},{3.0,2.0,3.0,12.0,-1.0},{4.0,-3.0,-5.0,-1.0,15.0} };
	z0 = { 1.0,1.0,1.0,1.0,1.0 };
	cout << endl << "Gauss-Seidel������" << endl;
	Gauss_Seidel_iteration(dim, A, z0, z, error, times);
	cout << "��������Ϊ = " << times << endl;
	cout << "������Ϊ��";
	Vector_cout(z0);

	A = { {10.0,1.0,2.0,3.0,4.0},{1.0,9.0,-1.0,2.0,-3.0},{2.0,-1.0,7.0,3.0,-5.0},{3.0,2.0,3.0,12.0,-1.0},{4.0,-3.0,-5.0,-1.0,15.0} };
	z0 = { 1.0,1.0,1.0,1.0,1.0 };
	cout << endl << "�����ݶȷ���" << endl;
	Conjugate_Gradient(dim, A, z, epsilon, z0, times, iteration_error);
	cout << "��������Ϊ = " << times << endl;
	cout << "�ο����Ϊ = " << iteration_error << endl;
	cout << "������Ϊ��";
	Vector_cout(z0);
	
}