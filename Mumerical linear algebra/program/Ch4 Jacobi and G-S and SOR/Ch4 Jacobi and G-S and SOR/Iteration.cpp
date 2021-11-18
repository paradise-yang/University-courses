#include "Iteration.h"

using namespace std;

typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

int main()
{
	//Ch4.1��Jacobi��Gauss-Seidel��SOR�����ɢ��ʽ΢�ַ��̣�epsilon=1.0 / 0.1 / 0.01 / 0.0001��
	
	double epsilon = 1.0;     //�˴��޸�epsilonֵ��
	double w = 1.75;           //SOR�������ɳ����ӣ�
	double error = 1e-5;
	clock_t start_time, end_time;
	int dim1 = 100, times = 0;
	
	double a = 0.5, h = 1.0 / dim1;//����4λ��Ч���֣��۲쾫ȷ�⣬����ȡ
	double ah2 = a * h * h;
	Matrix matrix(dim1 - 1, Vector(dim1 - 1, 0.0)), matrix0;
	Vector x0(dim1 - 1, 1.0), value(dim1 - 1, ah2), exact(dim1), y(dim1 - 1, 0);
	value[dim1 - 2] = value[dim1 - 2] - epsilon - h;
	
	cout << endl << "Ch4.1��Jacobi��Gauss-Seidel��SOR�����ɢ��ʽ΢�ַ��̣�epsilon = " << epsilon << endl;
	F(dim1, epsilon, a, y);   //��⾫ȷ�⣻
	cout << "��ȷ��Ϊ��" << endl;
	Vector_cout(y);

	cout << endl << "Jacobi������";
	band_matrix(dim1 - 1, matrix, epsilon, h);
	matrix0 = matrix;
	exact = x0;
	start_time = clock();                //��ȡ��һ��ʱ�䣻
	Jacobi_iteration(dim1 - 1, matrix, x0, value, error, times);
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "������Ϊ��" << endl;
	Vector_cout(x0);
	cout << "��������Ϊ = " << times << endl;
	cout << "����ʱ��Ϊ =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	Vector_subtraction(x0, y, x0);
	cout << "�뾫ȷ�����||x - x*|| = " << vector_infinity_norm(x0) << endl;
	cout << "�뾫ȷ�����||x - x*||2 = " << vector_2_norm(x0) << endl;
	matrix.clear();
	x0.clear();

	cout << endl << "Gauss-Seidel������";
	matrix = matrix0;
	x0 = exact;
	start_time = clock();                //��ȡ��һ��ʱ�䣻
	Gauss_Seidel_iteration(dim1 - 1, matrix, x0, value, error, times);
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "������Ϊ��" << endl;
	Vector_cout(x0);
	cout << "��������Ϊ = " << times << endl;
	cout << "����ʱ��Ϊ =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	Vector_subtraction(x0, y, x0);
	cout << "�뾫ȷ�����||x - x*|| = " << vector_infinity_norm(x0) << endl;
	cout << "�뾫ȷ�����||x - x*||2 = " << vector_2_norm(x0) << endl;
	matrix.clear();
	x0.clear();

	cout << endl << "SOR������" << endl << "�ɳ�����Ϊ w = " << w << endl;
	matrix = matrix0;
	x0 = exact;
	start_time = clock();                //��ȡ��һ��ʱ�䣻
	SOR_iteration(dim1 - 1, w, matrix, x0, value, error, times);
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "������Ϊ��" << endl;
	Vector_cout(x0);
	cout << "��������Ϊ = " << times << endl;
	cout << "����ʱ��Ϊ =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	Vector_subtraction(x0, y, x0);
	cout << "�뾫ȷ�����||x - x*|| = " << vector_infinity_norm(x0) << endl;
	cout << "�뾫ȷ�����||x - x*||2 = " << vector_2_norm(x0) << endl;
	matrix.clear();
	x0.clear();
	
	//Ch4.2 Jacobi��Gauss-Seidel��SOR�������ƫ΢�ַ��̣�N=20 / 40 /80��
	error = 1e-7;
	cout << endl << "Ch4.2 Guass-Seidel�������ƫ΢�ַ��̣�N=20 / 40 /80��" << endl;
	
	int N = 20;
	cout << "N = 20 ��" << endl;
	
	start_time = clock();                //��ȡ��һ��ʱ�䣻
	Jacobi_PDE(N, times, error);
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "Jacobi�������� = " << times;
	cout << "           ����ʱ��Ϊ =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	
	start_time = clock();                //��ȡ��һ��ʱ�䣻
	Gauss_Seidel_PDE(N, times, error);
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "Gauss-Seidel�������� = " << times;
	cout << "     ����ʱ��Ϊ =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;

	start_time = clock();                //��ȡ��һ��ʱ�䣻
	SOR_PDE(N, times, error, w);
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "SOR�������� = " << times;
	cout << "               ����ʱ��Ϊ =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl << endl;
	
	N = 40;
	cout << "N = 40 ��" << endl;

	start_time = clock();                //��ȡ��һ��ʱ�䣻
	Jacobi_PDE(N, times, error);
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "Jacobi�������� = " << times;
	cout << "          ����ʱ��Ϊ =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;

	start_time = clock();                //��ȡ��һ��ʱ�䣻
	Gauss_Seidel_PDE(N, times, error);
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "Gauss-Seidel�������� = " << times;
	cout << "     ����ʱ��Ϊ =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;

	start_time = clock();                //��ȡ��һ��ʱ�䣻
	SOR_PDE(N, times, error, w);
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "SOR�������� = " << times;
	cout << "              ����ʱ��Ϊ =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl << endl;
	
	N = 80;
	cout << "N = 80 ��" << endl;

	start_time = clock();                //��ȡ��һ��ʱ�䣻
	Jacobi_PDE(N, times, error);
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "Jacobi�������� = " << times;
	cout << "          ����ʱ��Ϊ =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;

	start_time = clock();                //��ȡ��һ��ʱ�䣻
	Gauss_Seidel_PDE(N, times, error);
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "Gauss-Seidel�������� = " << times;
	cout << "    ����ʱ��Ϊ =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;

	start_time = clock();                //��ȡ��һ��ʱ�䣻
	SOR_PDE(N, times, error, w);
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "SOR�������� = " << times;
	cout << "              ����ʱ��Ϊ =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl << endl;
	
}