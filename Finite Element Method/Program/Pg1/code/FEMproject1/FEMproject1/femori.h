#pragma once

#include <iomanip>
#include<vector>
#include<cmath>

using namespace std;

//����ֲ���һ�α�׼�����ڻ����þֲ���׼��������
void standard_product_1(vector<vector<double>>& a)
{
	a[0][0] = 1;
	a[0][1] = -1;
	a[1][1] = 1;
	a[1][0] = a[0][1];
}

//����ֲ��Ķ��α�׼�����ڻ����þֲ���׼��������
void standard_product_2(vector<vector<double>>& a)
{
	a[0][0] = 7.0 / 3, a[0][1] = -8.0 / 3, a[0][2] = 1.0 / 3;
	a[1][0] = a[0][1], a[1][1] = 16.0 / 3, a[1][2] = -8.0 / 3;
	a[2][0] = a[0][2], a[2][1] = a[1][2], a[2][2] = 7.0 / 3;
}

//�ֲ���ָ���Ӧ�����ָ�ꣻ
int index(int e, int sign)
{
	if (sign == 0) return e - 2;
	if (sign == 1) return e - 1;
}

//����һ�κ���������ϵ������K;
void make_cof_k(int N, vector<vector<double>>& cof_k, vector<vector<double>> a)
{
	//С���������N����ϵ������cof_k��СΪ(N-1)*(N-1)��
	double he = 1.0 / N;//��ɢ���䳤�ȣ�
	cof_k[0][0] = a[1][1] / he;
	for (int e = 2; e < N; e++)
	{
		cof_k[index(e, 0)][index(e, 0)] += (a[0][0] / he);
		cof_k[index(e, 1)][index(e, 1)] += (a[1][1] / he);
		cof_k[index(e, 0)][index(e, 1)] += (a[0][1] / he);
		cof_k[index(e, 1)][index(e, 0)] += (a[1][0] / he);
	}
	cof_k[N - 2][N - 2] += a[0][0] / he;
}

//���εľֲ���ָ���Ӧ�����ָ�ꣻ
int index_2(int e, int sign)
{
	if (sign == 0) return e - 2;
	if (sign == 1) return e - 1;
	if (sign == 2) return e;
}


//���ڶ��κ���������ϵ������K��
void make_cof_k_2(int N, vector<vector<double>>& cof_k, vector<vector<double>> a)
{
	//С���������N����ϵ������cof_k��СΪ(2N-1)*(2N-1)��
	double he = 1.0 / N;//��ɢ���䳤�ȣ�
	cof_k[0][0] = a[1][1] / he;
	cof_k[0][1] = a[1][2] / he;
	cof_k[1][0] = a[2][1] / he;
	cof_k[1][1] = a[2][2] / he;
	int e = 3;
	while (e < 2 * N - 2)
	{
		cof_k[index_2(e, 0)][index_2(e, 0)] += (a[0][0] / he);
		cof_k[index_2(e, 0)][index_2(e, 1)] += (a[0][1] / he);
		cof_k[index_2(e, 0)][index_2(e, 2)] += (a[0][2] / he);
		cof_k[index_2(e, 1)][index_2(e, 0)] += (a[1][0] / he);
		cof_k[index_2(e, 1)][index_2(e, 1)] += (a[1][1] / he);
		cof_k[index_2(e, 1)][index_2(e, 2)] += (a[1][2] / he);
		cof_k[index_2(e, 2)][index_2(e, 0)] += (a[2][0] / he);
		cof_k[index_2(e, 2)][index_2(e, 1)] += (a[2][1] / he);
		cof_k[index_2(e, 2)][index_2(e, 2)] += (a[2][2] / he);
		e += 2;
	}
	/*for (int e = 2; e < N; e++)
	{
		cof_k[index(e, 0)][index(e, 0)] += (a[0][0] / he);
		cof_k[index(e, 0)][index(e, 1)] += (a[0][1] / he);
		cof_k[index(e, 0)][index(e, 2)] += (a[0][2] / he);
		cof_k[index(e, 1)][index(e, 0)] += (a[1][0] / he);
		cof_k[index(e, 1)][index(e, 1)] += (a[1][1] / he);
		cof_k[index(e, 1)][index(e, 2)] += (a[1][2] / he);
		cof_k[index(e, 2)][index(e, 0)] += (a[2][0] / he);
		cof_k[index(e, 2)][index(e, 1)] += (a[2][1] / he);
		cof_k[index(e, 2)][index(e, 2)] += (a[2][2] / he);
	}*/
	cof_k[2 * N - 3][2 * N - 3] += a[0][0] / he;
	cof_k[2 * N - 3][2 * N - 2] += a[0][1] / he;
	cof_k[2 * N - 2][2 * N - 3] += a[1][0] / he;
	cof_k[2 * N - 2][2 * N - 2] += a[1][1] / he;
}

//����f(x)��
double f(double x)
{
	return -2 * cos(x) + (x - 1) * sin(x);
}

//��i��һ�λ�������
double basis_1(int N, int i, double x)
{
	double he = 1.0 / N;
	if ((x > (i - 1.0) * he) && (x <= i * he))
		return (x - (i - 1.0) * he) / he;
	else if ((x > i * he) && (x < (i + 1.0) * he))
		return ((i + 1.0) * he - x) / he;
	else
		return 0.0;
}

//��i�����λ�������
double basis_2(int N, int i, double x)
{
	double he = 1.0 / N;
	switch (i % 2)
	{
	case 0:
		if ((x > (i / 2 - 1) * he) && (x <= i / 2 * he))
			return (2 * x - (i - 1) * he) * (x - (i / 2 - 1) * he) / he / he;
		else if ((x > i / 2 * he) && (x < (i / 2 + 1) * he))
			return (2 * x - (i + 1) * he) * (x - (i / 2 + 1) * he) / he / he;
		else
			return 0;
	case 1:
		if ((x > (i - 1) / 2 * he) && x < (i + 1) / 2 * he)
			return 4 * (x - (i - 1) / 2 * he) * ((i + 1) / 2 * he - x) / he / he;
		else
			return 0;
	}
	/*
	if (i % 2 ==1)
	{
		if ((x > (i - 1) * he) && x < (i + 1) * he)
			return -1.0 / he / he * pow(x - i * he, 2) + 1;
		else
			return 0;
	}
	else
	{
		if ((x > (i / 2 - 1) * he) && (x <= i / 2 * he))
			return (2 * x - (i - 1) * he) * (x - (i / 2 - 1) * he) / he / he;
		else if ((x > i/2 * he) && (x < (i/2+1) * he))
			return (1.0 / he * (x - i * he) - 1) * (1.0 / 2 / he * (x - i * he) - 1);
		else
			return 0;
	}
	*/
}

//����һ�κ����ĸ�������Gauss���ֹ�ʽ��
double gauss_integrate_1(int N, int i, double left, double right)
{
	int n = 16;
	double delta_x = (right - left) / n;
	double temp = 0;
	for (int j = 0; j < n; j++)
	{
		double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
		double x1 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 - sqrt(15) / 10 * delta_x / 2;
		double x2 = (left + j * delta_x + left + (j + 1) * delta_x) / 2;
		double x3 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 + sqrt(15) / 10 * delta_x / 2;
		temp+= a1 * f(x1) * basis_1(N, i, x1) + a2 * f(x2) * basis_1(N, i, x2) + a3 * f(x3) * basis_1(N, i, x3);
	}
	return temp;
}

//���ڶ��κ����ĸ�������Gauss���ֹ�ʽ��
double gauss_integrate_2(int N, int i, double left, double right)
{
	int n = 16;
	double delta_x = (right - left) / n;
	double temp = 0;
	for (int j = 0; j < n; j++)
	{
		double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
		double x1 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 - sqrt(15) / 10 * delta_x / 2;
		double x2 = (left + j * delta_x + left + (j + 1) * delta_x) / 2;
		double x3 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 + sqrt(15) / 10 * delta_x / 2;
		temp += a1 * f(x1) * basis_2(N, i, x1) + a2 * f(x2) * basis_2(N, i, x2) + a3 * f(x3) * basis_2(N, i, x3);
	}
	return temp;
}

//�����Ҳ�������
void make_cof_f(int N, vector<double>& cof_f_1, vector<double>& cof_f_2)
{
	double he = 1.0 / N;
	for (int i = 1; i < N; i++)
	{
		cof_f_1[i - 1] = gauss_integrate_1(N, i, (i - 1.0) * he, (i + 1.0) * he);
		if (i % 2 == 0)
			cof_f_2[i - 1] = gauss_integrate_2(N, i, (i / 2 - 1) * he, (i / 2 + 1) * he);
		else
			cof_f_2[i - 1] = gauss_integrate_2(N, i, (i - 1) / 2 * he, (i + 1) / 2 * he);
	}
	for (int i = N; i < 2 * N; i++)
	{
		if (i % 2 == 0)
			cof_f_2[i - 1] = gauss_integrate_2(N, i, (i / 2 - 1) * he, (i / 2 + 1) * he);
		else
			cof_f_2[i - 1] = gauss_integrate_2(N, i, (i - 1) / 2 * he, (i + 1) / 2 * he);
	}
}

//���������������
void Matrix_times_vector(vector<vector<double>> matrix, vector<double> column, vector<double>& value)
{
	int dim = size(matrix);
	for (int i = 0; i < dim; i++)
	{
		value[i] = 0;
		for (int j = 0; j < dim; j++)
			value[i] += matrix[i][j] * column[j];
	}
}

//����֮����������
void Vector_subtraction(vector<double> x, vector<double> y, vector<double>& dx)
{
	int dim = 0;
	dim = size(x);
	for (int i = 0; i < dim; i++)
		dx[i] = x[i] - y[i];
}

//�������ڻ���
double Inner_product_of_vectors(vector<double> x, vector<double> y)
{
	int i, n;
	double temp = 0;
	n = size(x);
	for (i = 0; i < n; i++)
		temp += x[i] * y[i];
	return temp;
}

//������2-�����������ط���ֵ��
double vector_2_norm(vector<double> x)
{
	int dim = size(x);
	double temp = 0.0;
	for (int i = 0; i < dim; i++)
		temp = temp + x[i] * x[i];
	return sqrt(temp);
}

//��������ptAp��
double pAp(vector<vector<double>> matrix, vector<double> column)
{
	int dim = size(matrix);
	vector<double> value(dim, 0);
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			value[i] += matrix[i][j] * column[j];
	return Inner_product_of_vectors(column, value);
}

//���������������ݶȵ��������µ�xk��
void xk(vector<double>& x, int dim, vector<double> x0, double a0, vector<double> p0)
{
	for (int i = 0; i < dim; i++)
		x[i] = x0[i] + a0 * p0[i];
}

//�������������ɲ�����rk��
void rk(vector<double>& r, int dim, vector<double> r0, vector<double> p0, double a0, vector<vector<double>> matrix)
{
	vector<double> temp(dim, 0);
	Matrix_times_vector(matrix, p0, temp);
	for (int i = 0; i < dim; i++)
		r[i] = r0[i] - a0 * temp[i];
}

//�����ݶȷ���
void Conjugate_Gradient(int dim, vector<vector<double>> matrix, vector<double> b, double epsilon, vector<double>& x, int& times, double& error)
{
	double e = epsilon * vector_2_norm(b);
	vector<double> x0(dim, 1);
	vector<double> r0(dim, 0), r(dim, 0), p0(dim, 0);
	Matrix_times_vector(matrix, x0, r0);
	Vector_subtraction(b, r0, r0);//˳��
	r = r0;
	x = x0;
	times = 0;
	while (vector_2_norm(r) > e && times < 10000)
	{
		times++;
		if (times == 1) p0 = r0;
		else
		{
			double beta = Inner_product_of_vectors(r, r) / Inner_product_of_vectors(r0, r0);
			xk(p0, dim, r, beta, p0);
		}
		double alpha = Inner_product_of_vectors(r, r) / pAp(matrix, p0);
		xk(x, dim, x, alpha, p0);
		r0 = r;//����rk-1;
		rk(r, dim, r0, p0, alpha, matrix);
	}
	error = vector_2_norm(r);
}

//���Խ���Gauss��ȥ��ⷽ���飻
void gauss(vector<vector<double>>& cof_k, vector<double>& cof_f, vector<double>& cof_u)
{
	//��ȥ�����ǲ��޸��Ҳ�ϵ����
	for (int i = 1;i<cof_k.size(); i++)
	{
		cof_k[i][i] = cof_k[i][i] - cof_k[i][i - 1] / cof_k[i - 1][i - 1] * cof_k[i - 1][i];
		cof_f[i] = cof_f[i] - cof_k[i][i - 1] / cof_k[i - 1][i - 1] * cof_f[i - 1];
		cof_k[i][i - 1] = 0;
	}

	//��ⷽ���飻
	for (int i = cof_k.size() - 1; i > 0; i--)
	{
		cof_u[i] = cof_f[i] / cof_k[i][i];
		//cout << cof_f[i] <<"  "<< cof_k[i][i] <<"  "<< cof_u[i] << endl;
		cof_f[i - 1] -= cof_u[i] * cof_k[i - 1][i];
	}
	cof_u[0] = cof_f[0] / cof_k[0][0];
}

//��ʵ�⺯��u(x)��
double u(double x)
{
	return (x - 1) * sin(x);
}

//����1�λ���������⺯��uh(x)��
double uh_basis_1(int N, vector<double> cof_u, double x)
{
	double temp = 0;
	double he = 1.0 / N;
	for (int i = 1; i < N; i++)
		temp += cof_u[i - 1] * basis_1(N, i, x);
	return temp;
}

//��ɢ�����������
double error_max_1(int N, vector<double> cof_u)
{
	int n = 100;
	double delta_x = 1.0 / n;
	vector<double> error(n + 1, 0);
	for (int i = 0; i <= n; i++)
		error[i] = fabs(u(i * delta_x) - uh_basis_1(N, cof_u, i * delta_x));
	double max = error[0];
	for (int i = 1; i <= n; i++)
		if (error[i] > max)
			max = error[i];
	return max;
}


//���ڸ���Gauss���ֹ�ʽ��ɢ����L1������
double error_L1_1(int N, vector<double> cof_u)
{
	int n = 128;
	double delta_x = 1.0 / n;
	double temp = 0;
	for (int j = 0; j < n; j++)
	{
		double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
		double x1 = (j * delta_x + (j + 1) * delta_x) / 2 - sqrt(15) / 10 * delta_x / 2;
		double x2 = (j * delta_x + (j + 1) * delta_x) / 2;
		double x3 = (j * delta_x + (j + 1) * delta_x) / 2 + sqrt(15) / 10 * delta_x / 2;
		temp += a1 * fabs(u(x1) - uh_basis_1(N, cof_u, x1)) + a2 * fabs(u(x2) - uh_basis_1(N, cof_u, x2)) + a3 * fabs(u(x3) - uh_basis_1(N, cof_u, x3));
	}
	return temp;
}

//����2�λ���������⺯��uh(x)��
double uh_basis_2(int N, vector<double> cof_u, double x)
{
	double temp = 0.0;
	double he = 1.0 / N;
	for (int i = 1; i < (2 * N); i++)
		temp += cof_u[i - 1] * basis_2(N, i, x);
	return temp;
}

//��ɢ�����������
double error_max_2(int N, vector<double> cof_u)
{
	int n = 100;
	double delta_x = 1.0 / n;
	vector<double> error(n + 1, 0);
	//cout << endl << "N=" << N << "��ֵ��" << endl;
	for (int i = 0; i <= n; i++)
		//cout << "{" << i * delta_x << "," << uh_basis_2(N, cof_u, i * delta_x) << "},",
		error[i] = fabs(u(i * delta_x) - uh_basis_2(N, cof_u, i * delta_x));
	double max = error[0];
	for (int i = 1; i <= n; i++)
		if (error[i] > max)
			max = error[i];
	return max;
}

//������������⣬����������
void error_out(vector<vector<double>>& error)
{
	//�����׼�����ڻ���
	vector<vector<double>> a1(2, vector<double>(2, 0)), a2(3, vector<double>(3, 0));
	standard_product_1(a1);
	standard_product_2(a2);

	for (int i = 0; i < 5; i++)
	{
		//������ɢ�̶ȣ�
		int N = pow(2, i) * 10;
		double he = 1.0 / N;
		//ϵ������K��
		vector<vector<double>> cof_k_1(N - 1, vector<double>(N - 1, 0)), cof_k_2(2 * N - 1, vector<double>(2 * N - 1, 0));
		make_cof_k(N, cof_k_1, a1);
		make_cof_k_2(N, cof_k_2, a2);
		//�Ҳຯ��f����ڻ�������
		vector<double> cof_f_1(N - 1, 0), cof_f_2(2 * N - 1, 0);
		make_cof_f(N, cof_f_1, cof_f_2);
		//��ɢ���ϵ��������
		vector<double> cof_u_1(N - 1, 0), cof_u_2(2 * N - 1, 1);
		///*
		//������
		if (N == 10)
		{
			for (int i = 0; i < size(cof_k_2); i++)
			{
				for (int j = 0; j < cof_k_2.size(); j++)
					cout << cof_k_2[i][j] << "  ";
				cout << endl;
			}
			cout << endl << endl;
			for (int i = 0; i < cof_f_2.size(); i++)
				cout << cof_f_2[i] << "  ";
			cout << endl << endl;
		}
		//*/
		//�����ݶȷ���⣻
		int times = 0;             //�����ݶȵ���������������
		double iteration_error = 0;//�����ݶȵ�����������ο���
		gauss(cof_k_1, cof_f_1, cof_u_1);
		//Conjugate_Gradient(N - 1, cof_k_1, cof_f_1, 1e-10, cof_u_1, times, iteration_error);
		times = 0;           //�����ݶȵ���������������
		iteration_error = 0;//�����ݶȵ�����������ο���
		Conjugate_Gradient(2 * N - 1, cof_k_2, cof_f_2, 1e-10, cof_u_2, times, iteration_error);
		//�����������Ϊ��
		error[i][0] = error_max_1(N, cof_u_1);
		error[i][2] = error_L1_1(N, cof_u_1);
		if (i > 0)
		{
			error[i][1] = log(error[i - 1][0] / error[i][0]) / log(2);
			error[i][3] = log(error[i - 1][2] / error[i][2]) / log(2);
		}
	}
}