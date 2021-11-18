#include <iostream>
#include<iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <ctime>

using namespace std;

//���������
void Matrix_cout(vector<vector<double>> matrix, int m, int n)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			cout << setw(9) << matrix[i][j] << "    ";
		cout << endl;
	}
}

//����˷�house*matrix,�������ڵ������������棻
void Matrix_times(vector<vector<double>> house, vector<vector<double>> matrix, vector<vector<double>>& A)
{
	int dim = size(matrix);
	vector<vector<double>> temp(dim, vector<double>(dim, 0));
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			for (int k = 0; k < dim; k++)
				temp[i][j] += house[i][k] * matrix[k][j];
		}
	}
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			A[i][j] = temp[i][j];
}

//���ɶ���ʽ���������棩�ѷ���
void Companion_matrix(vector<vector<double>>& matrix, vector<double> poly)
{
	int dim = size(poly);
	int n = dim - 1;
	for (int i = 0; i < n; i++)
		matrix[i][n - 1] = 0 - poly[i];
	for (int i = 1; i < n; i++)
		matrix[i][i - 1] = 1.0;
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

//������������������ط���ֵ��
double vector_infinity_norm(vector<double> x)
{
	int i, n;
	double norm = fabs(x[0]);
	n = size(x);
	for (i = 1; i < n; i++)
	{
		if (fabs(x[i]) > norm)
			norm = fabs(x[i]);
	}
	return norm;
}

//������ģ��������
double magnitude_max_component_of_vector(vector<double> x)
{
	int dim = size(x);
	double norm = fabs(x[0]);
	int temp = 0;
	for (int i = 1; i < dim; i++)
	{
		if (fabs(x[i]) > norm)
		{
			norm = fabs(x[i]);
			temp = i;
		}
	}
	return x[temp];
}

//����֮����������
void Vector_subtraction(vector<double> x, vector<double> y, vector<double>& dx)
{
	int dim = 0;
	dim = size(x);
	for (int i = 0; i < dim; i++)
		dx[i] = x[i] - y[i];
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

//������������yk��muk����uk��
void uk(vector<double> yk, double muk, vector<double>& uk)
{
	int dim = size(yk);
	for (int i = 0; i < dim; i++)
		uk[i] = yk[i] / muk;
}

//�ݷ���ⷽ��ģ�������ֵ��
void Power_method(int dim, vector<vector<double>> A, double error , double& value)
{
	vector<double> u0(dim, 1);//��ʼ������
	vector<double> y(dim, 0), u(dim, 1), du(dim, 1);
	double mu = 0;
	do {
		Matrix_times_vector(A, u0, y);//�������yk
		mu = magnitude_max_component_of_vector(y);//��ykģ��������
		uk(y, mu, u);
		Vector_subtraction(u, u0, du);
		u0 = u;
	} while (vector_infinity_norm(du) > error);
	//value = mu;//���õ����ó���mu����ģ�������ֵ���ձ����ã�
	value = A[0][dim - 1] * u[dim - 1] / u[0];//���������������ѷ������ԣ��Լ�����ֵ���壬���õ�һ�������������ֵ��ֻ�����ѷ������Σ�
}

//�ݷ������ڶ���ʽ�ѷ������ģ������
double max_modulus_root_polynomial(vector<double> poly, double error)
{
	int dim = size(poly);
	int n = dim - 1;//����ʽ����������ά����
	vector<vector<double>> matrix(n, vector<double>(n, 0));
	//�����ѷ���
	for (int i = 0; i < n; i++)
		matrix[i][n - 1] = 0 - poly[dim - 1 - i];
	//for (int i = 0; i < n; i++)//��ʦ�Ͽν����ѷ���ʵ��ֻ�������һ������ֵ��������ģ���
	//	matrix[0][i] = 0 - poly[i + 1];
	for (int i = 1; i < n; i++)
		matrix[i][i - 1] = 1.0;
	//Matrix_cout(matrix, n, n);
	double value = 0;
	Power_method(n, matrix, error, value);
	return value;
}

//����Householder�任��ʹ��������Ϊe1��������������v��ϵ��b��
void Householder(vector<double> x, vector<double>& v, double& b)
{
	int dim = size(x);
	double max = 0, m = 0, a = 0;
	max = vector_infinity_norm(x);
	x[0] = x[0] / max;
	for (int i = 1; i < dim; i++)
	{
		v[i] = x[i] / max;
		m += pow(v[i], 2);
	}
	if (m == 0) b = 0;
	else
	{
		a = sqrt(x[0] * x[0] + m);
		if (x[0] <= 0) v[0] = x[0] - a;
		else v[0] = -m / (x[0] + a);
		b = 2 * pow(v[0], 2) / (m + pow(v[0], 2));
		for (int i = 1; i < dim; i++)
			v[i] = v[i] / v[0];
		v[0] = 1.0;
	}
}

//Householder�������ɣ�
void Householder_make(int n, vector<vector<double>>& house, vector<double> v, double b)
{
	int dim = size(house);
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			if (i > n && j > n) 
			{
				if (i == j) house[i][j] = 1.0 - b * v[i - n - 1] * v[j - n - 1];
				else house[i][j] = 0 - b * v[i - n - 1] * v[j - n - 1];
			}
			else
			{
				if (i == j) house[i][j] = 1;
				else house[i][j] = 0;
			}
		}
	}
}

//6.4.1����������Householder���������ã�
void householder_opertion(vector<vector<double>>& A, int n, vector<double> v,double b)
{
	int dim = size(A);
	int vl = size(v);
	vector<vector<double>> temp(vl, vector<double>(vl, 0));
	for (int i = 0; i < vl; i++)
		for (int j = 0; j < vl; j++)
			for (int k = 0; k < vl; k++)
				temp[i][j] += b * v[i] * v[j] * A[n + 1 + k][n + 1 + j];
	for (int i = 0; i < vl; i++)
		for (int j = 0; j < vl; j++)
			A[n + 1 + i][n + 1 + j] -= temp[i][j];
	for (int i = 2; i < dim; i++)
		for (int j = 0; j < i - 1; j++)
			A[i][j] = 0;

	vector<vector<double>> temp0(vl, vector<double>(vl, 0));
	for (int i = 0; i < vl; i++)
		for (int j = 0; j < vl; j++)
			for (int k = 0; k < vl; k++)
				temp0[i][j] += A[n + 1 + i][n + 1 + k] * b * v[k] * v[j];
	for (int i = 0; i < vl; i++)
		for (int j = 0; j < vl; j++)
			A[n + 1 + i][n + 1 + j] -= temp0[i][j];
	for (int i = 2; i < dim; i++)
		for (int j = 0; j < i - 1; j++)
			A[i][j] = 0;
}

//�㷨6.4.1 ��Householder��������������Hessenberg�ֽ⣬�ֽⴢ����A������任���󴢴���Q�
void Hessenberg_decomposition(vector<vector<double>>& A)
{
	int dim = size(A);
	for (int j = 0; j < dim - 2; j++)
	{
		double b = 0;
		int tempk = dim - j - 1;
		vector<double> v(tempk, 0);
		vector<double> x(tempk, 0);  //��ȡi�жԽ�Ԫ��������Householder�任��
		for (int i = 0; i < tempk; i++)
			x[i] = A[i + j + 1][j];
		Householder(x, v, b);        //��Householder�任��
		vector<vector<double>> house(dim, vector<double>(dim, 0));
		Householder_make(j, house, v, b);//����Householder�任��
		Matrix_times(house, A, A);//Householder�任��ˣ�
		for (int i = j + 2; i < dim; i++) A[i][j] = 0;
		Matrix_times(A, house, A);//Householder�任�ҳˣ�
		
	}
	//������������
	for (int i = 2; i < dim; i++)
		for (int j = 0; j < i - 1; j++)
			A[i][j] = 0;
}

//����Hesssenberg����СԪ����Ϊ0��
void clean_small_elements(int dim, vector<vector<double>>& H, double u)
{
	for (int i = 1; i < dim; i++)
		if (fabs(H[i][i - 1]) <= (fabs(H[i][i]) + fabs(H[i - 1][i - 1])) * u) H[i][i - 1] = 0;
}

//ȷ����ʽQR���ᵽ��m,l������ע��mΪ���е�m������������Ĭ��H����Hessenberg��
void confirm_m_l(vector<vector<double>> H, int& m, int& l)
{
	int n = H.size(), t = 1;
	while ((m < n - 2) && (t == 1))
	{
		if (H[n - m - 1][n - m - 2] == 0) m++;
		else if ((H[n - m - 2][n - m - 3] == 0) && ((H[n - m - 2][n - m - 2] - H[n - m - 1][n - m - 1]) * (H[n - m - 2][n - m - 2] - H[n - m - 1][n - m - 1]) + 4 * H[n - m - 2][n - m - 1] * H[n - m - 1][n - m - 2]) < 0)
			//�ж�2������
		{
			m += 2, t += 2;
		}
		else t = 0;
	}
	if (m == n - 1) m = n, l = 0;
	else if (m == n - 2)
	{
		if (((H[0][0] - H[1][1]) * (H[0][0] - H[1][1]) + 4 * H[0][1] * H[1][0] < 0)) m = n, l = 0;
		else l = 2;
	}

	if (m != n) for (l = 2; (H[n - m - l + 1][n - m - l] != 0) && (l < (n - m)); l++);
	l = n - m - l;
}

//˫�ز�λ�Ƶ�QR����&Ĭ��PΪ0��
void Francis(vector<vector<double>>& H, vector<vector<double>>& P)
{
	int n = H.size(), q, r, i;
	double s, t, x, y, z, a, b, c, d;
	vector<double> v(3), w(3);
	for (i = 0; i < n; i++) P[i][i] = 1;
	if (n >= 3)
	{
		s = H[n - 2][n - 2] + H[n - 1][n - 1];
		t = H[n - 2][n - 2] * H[n - 1][n - 1] - H[n - 2][n - 1] * H[n - 1][n - 2];
		x = H[0][0] * H[0][0] + H[0][1] * H[1][0] - s * H[0][0] + t;
		y = H[1][0] * (H[0][0] + H[1][1] - s);
		z = H[1][0] * H[2][1];
		for (int k = 0; k < n - 2; k++)
		{
			a = sqrt(x * x + y * y + z * z);
			v[0] = x - a, v[1] = y, v[2] = z;
			b = a * a - a * x;                       //��Householder�任�е�beta���ʽ��
			if (b != 0)
			{
				b = 1 / b;
				q = (1 > k) ? 1 : k;
				w[0] = b * v[0], w[1] = b * v[1], w[2] = b * v[2];
				for (i = q - 1; i < n; i++)
				{
					c = v[0] * H[k][i] + v[1] * H[k + 1][i] + v[2] * H[k + 2][i];
					H[k][i] -= w[0] * c, H[k + 1][i] -= w[1] * c, H[k + 2][i] -= w[2] * c;
				}
				r = (k + 4 > n) ? n : k + 4;
				for (i = 0; i < r; i++)
				{
					c = H[i][k] * v[0] + H[i][k + 1] * v[1] + H[i][k + 2] * v[2];
					H[i][k] -= c * w[0], H[i][k + 1] -= c * w[1], H[i][k + 2] -= c * w[2];
					d = P[i][k] * v[0] + P[i][k + 1] * v[1] + P[i][k + 2] * v[2];
					P[i][k] -= d * w[0], P[i][k + 1] -= d * w[1], P[i][k + 2] -= d * w[2];
				}
				x = H[k + 1][k], y = H[k + 2][k];
				if (k < n - 3) z = H[k + 3][k];
			}
		}
		a = sqrt(x * x + y * y);
		v[0] = x - a, v[1] = y;
		b = a * a - a * x;
		if (b != 0)
		{
			b = 1 / b;
			w[0] = b * v[0], w[1] = b * v[1];
			for (i = n - 3; i < n; i++)
			{
				c = v[0] * H[n - 2][i] + v[1] * H[n - 1][i];
				H[n - 2][i] -= w[0] * c, H[n - 1][i] -= w[1] * c;
			}
			for (i = 0; i < n; i++)
			{
				c = H[i][n - 2] * v[0] + H[i][n - 1] * v[1];
				H[i][n - 2] -= w[0] * c, H[i][n - 1] -= w[1] * c;
				d = P[i][n - 2] * v[0] + P[i][n - 1] * v[1];
				P[i][n - 2] -= w[0] * d, P[i][n - 1] -= w[1] * d;
			}
		}
	}
	else if (n == 2)
		//��ͨ��Judgegment���������ķ����Ȼ����s�Ǹ�������P������Ϊ��λ��
	{
		a = H[0][0], b = H[0][1], c = H[1][0], d = H[1][1];
		s = (a - d) * (a - d) + 4 * b * c;
		if (s >= 0)
		{
			H[0][1] = b - c, H[1][0] = 0, H[0][0] = (a + d - sqrt(s)) / 2, H[1][1] = (a + d + sqrt(s)) / 2;
		}
		if (b == 0) P[0][0] = 0, P[0][1] = 1, P[1][0] = -1, P[1][1] = 0;
		else t = (a - H[0][0]) / b, P[0][0] = P[1][1] = 1 / sqrt(1 + t * t), P[0][1] = P[0][0] * t, P[1][0] = -P[0][1];
	}
}

//6.4.3������������ȡ����
void part_make(vector<vector<double>> H, vector<vector<double>>& part, int starti, int endi, int startj, int endj)
{
	for (int i = starti; i <= endi; i++)
		for (int j = startj; j <= endj; j++)
			part[i - starti][j - startj] = H[i][j];
}

//6.4.3����������H12��H23��P�Ĳ�����
void submatrix_operations(vector<vector<double>>& H, vector<vector<double>> P, int m, int l)
{
	int n = H.size();
	vector<vector<double>> Q(l, vector<double>(n - m - l, 0)), R(n - m - l, vector<double>(m, 0));
	//H12�Ĳ���
	for (int i = 0; i < l; i++)
		for (int j = 0; j < n - m - l; j++)
			for (int k = 0; k < n - m - l; k++)
				Q[i][j] += H[i][k + l] * P[k][j];
	for (int i = 0; i < l; i++)
		for (int j = 0; j < n - m - l; j++)
			H[i][j + l] = Q[i][j];
	//H23�Ĳ���
	for (int i = 0; i < n - m - l; i++)
		for (int j = 0; j < m; j++)
			for (int k = 0; k < n - m - l; k++)
				R[i][j] += P[k][i] * H[l + k][n - m + j];
	for (int i = 0; i < n - m - l; i++)
		for (int j = 0; j < m; j++)
			H[l + i][n - m + j] = R[i][j];
}

//�㷨6.4.3 ����ʵ�����ʵSchur�ֽ⣬����ʽQR�㷨��
void Real_Schur_decomposition(vector<vector<double>>& A, double u)
{
	int dim = size(A);
	int m = 0, l = 0;
	Hessenberg_decomposition(A);        //Step2��
	clean_small_elements(dim, A, u);    //Step3(1)��
	confirm_m_l(A, m, l);
	while (m != dim)
	{
		//Step4��
		vector<vector<double>> H22(dim - m - l, vector<double>(dim - m - l, 0)), P(dim - m - l, vector<double>(dim - m - l, 0));
		part_make(A, H22, l, dim - m - 1, l, dim - m - 1);//ȡ��H22��
		Francis(H22, P);
		for (int i = l; i < dim - m; i++)                 //��ֵ��ȥ��
			for (int j = l; j < dim - m; j++)
				A[i][j] = H22[i - l][j - l];

		//Step5��
		submatrix_operations(A, P, m, l);
		for (int i = 1; i < dim; i++) 
			if (abs(A[i][i - 1]) < 0.00001) A[i][i - 1] = 0;

		//������������
		for (int i = 2; i < dim; i++)
			for (int j = 0; j < i - 1; j++)
				A[i][j] = 0;
		
		clean_small_elements(dim, A, u);//Step3(1)��
		confirm_m_l(A, m, l);           //Step3(2)��
	}
	for (int i = 1; i < dim; i++) if (abs(A[i][i - 1]) < 0.0000001) A[i][i - 1] = 0;
	
	/*��do while д�����ޱ��ʲ���,��������6.4.3��˳��
	do {
		clean_small_elements(dim, A, u);//Step3(1)��
		int m = 0, l = 0;
		confirm_m_l(A, m, l);           //Step3(2)��
		if (m >= dim) break;            //Step3(3)��

		//Step4��
		vector<vector<double>> H22(dim - m - l, vector<double>(dim - m - l, 0)), P(dim - m - l, vector<double>(dim - m - l, 0));
		part_make(A, H22, l, dim - m - 1, l, dim - m - 1);//ȡ��H22��
		Francis(H22, P);
		for (int i = l; i < dim - m; i++)                  //��ֵ��ȥ��
			for (int j = l; j < dim - m; j++)
				A[i][j] = H22[i - l][j - l];

		//Step5��
		submatrix_operations(A, P, m, l);
		for (int i = 1; i < dim; i++) if (abs(A[i][i - 1]) < 0.00001) A[i][i - 1] = 0;

		//������������
		for (int i = 2; i < dim; i++)
			for (int j = 0; j < i - 1; j++)
				A[i][j] = 0;
	} while (1);*/
}

//���������Schur�ֽ⣬��ⷽ�������ֵ�������
void Eigenvalue(vector<vector<double>> A)
{
	int dim = size(A);
	int i = 0;
	while (i < dim - 1)
	{
		if (A[i + 1][i] == 0)
			cout << setw(9) << A[i][i] << "  ", i++;
		else 
		{
			double a = A[i][i], b = A[i][i + 1], c = A[i + 1][i], d = A[i + 1][i + 1];
			cout << setw(9) << (a + d) / 2 << "��" << sqrt(-4 * b * c - (a - d) * (a - d)) / 2 << "i" << "  ", i += 2;
		}
	}
	if (i == dim - 1)
		cout << setw(9) << A[dim - 1][dim - 1] << "  ";
	cout << endl;
}