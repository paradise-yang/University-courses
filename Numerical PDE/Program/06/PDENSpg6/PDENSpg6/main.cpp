#include<iostream>
#include<iomanip>
#include<vector>
#include<cmath>
#include<fstream>
#include<string>

using namespace std;

#define pi 3.1415926535897932384626433832795028841971 //Բ����\piֵ��
#define miu 0.1 //����ȣ�

//��ȷ�⺯����
double u_exact(double x, double t)
{
	if (x >= -pi && x <= 0)
		return 0.5*exp(-4 * t)*sin(x);
	else
		return exp(-4 * t)*sin(2 * x);
}
//��ֵ������
double a(double x)
{
	if (x < 0)
		return 4.0;
	else
		return 1.0;
}
//FTBS��ʽ��
void FTCS_mean(int J, double lamuda, double T, vector<vector<double>> &u)
{
	double deltax = 2.0 * pi / J;              //�ռ���ɢ�̶ȣ�
	double deltat = lamuda * deltax * deltax;//���ռ���ɢ�̶ȣ�
	int N = floor(T / deltat);               //��Ҫǰ���Ĳ�����
	double last_deltat = T - N * deltat;
	//��ֵ��
	for (int j = 0; j <= J; j++)
		u[0][j] = u_exact(j * deltax - pi, 0);
	//FTBS��ʽ��
	for (int i = 1; i <= N; i++)
	{
		u[i % 2][0] = 0;//��߽紦��
		for (int j = 1; j < J; j++)
		{
			double aleft = (a((j - 1) * deltax - pi) + a(j * deltax - pi)) / 2;
			double aright= (a((j + 1) * deltax - pi) + a(j * deltax - pi)) / 2;
			u[i % 2][j] = lamuda * aleft*u[(i + 1) % 2][j - 1] + (1 - lamuda * aright - lamuda * aleft)*u[(i + 1) % 2][j] + lamuda * aright*u[(i + 1) % 2][j + 1];
		}
		u[i % 2][J] = 0;//�ұ߽紦��
	}
	//���һ���̲�������
	double lamuda_last = last_deltat / deltax / deltax;
	int i = N + 1;
	u[i % 2][0] = 0;//��߽紦��
	for (int j = 1; j < J; j++)
	{
		double aleft = (a((j - 1) * deltax - pi) + a(j * deltax - pi)) / 2;
		double aright = (a((j + 1) * deltax - pi) + a(j * deltax - pi)) / 2;
		u[i % 2][j] = lamuda * aleft*u[(i + 1) % 2][j - 1] + (1 - lamuda * aright - lamuda * aleft)*u[(i + 1) % 2][j] + lamuda * aright*u[(i + 1) % 2][j + 1];
	}
	u[i % 2][J] = 0;//�ұ߽紦��
}
//
double theta(double x, double xj_1, double xj, double deltax)
{
	if (x >= xj_1 && x <= xj)
		return 1.0 / 2;
		//return (x - xj) / deltax;
	else
		return (x - xj) / deltax;
		//return 1.0 / 2;
}
void FTCS_harmonic(int J, double lamuda, double T, vector<vector<double>> &u)
{
	double deltax = 2.0 * pi / J;              //�ռ���ɢ�̶ȣ�
	double deltat = lamuda * deltax * deltax;//���ռ���ɢ�̶ȣ�
	int N = floor(T / deltat);               //��Ҫǰ���Ĳ�����
	double last_deltat = T - N * deltat;
	//��ֵ��
	for (int j = 0; j <= J; j++)
		u[0][j] = u_exact(j * deltax - pi, 0);
	//FTBS��ʽ��
	for (int i = 1; i <= N; i++)
	{
		u[i % 2][0] = 0;//��߽紦��
		for (int j = 1; j < J; j++)
		{
			double aleft = 1.0 / (theta(0, (j - 1) * deltax - pi, j * deltax - pi, deltax) / a((j - 1) * deltax - pi) + (1.0 - theta(0, (j - 1) * deltax - pi, j * deltax - pi, deltax)) / a(j * deltax - pi));
				//(a((j - 1) * deltax - pi) + a(j * deltax - pi)) / 2;
			double aright = 1.0 / (theta(0, j * deltax - pi, (j + 1) * deltax - pi, deltax) / a(j * deltax - pi) + (1.0 - theta(0, j * deltax - pi, (j + 1) * deltax - pi, deltax)) / a((j + 1) * deltax - pi));
				//(a((j + 1) * deltax - pi) + a(j * deltax - pi)) / 2;
			u[i % 2][j] = lamuda * aleft*u[(i + 1) % 2][j - 1] + (1 - lamuda * aright - lamuda * aleft)*u[(i + 1) % 2][j] + lamuda * aright*u[(i + 1) % 2][j + 1];
		}
		u[i % 2][J] = 0;//�ұ߽紦��
	}
	//���һ���̲�������
	double lamuda_last = last_deltat / deltax / deltax;
	int i = N + 1;
	u[i % 2][0] = 0;//��߽紦��
	for (int j = 1; j < J; j++)
	{
		double aleft = 1.0 / (theta(0, (j - 1) * deltax - pi, j * deltax - pi, deltax) / a((j - 1) * deltax - pi) + (1.0 - theta(0, (j - 1) * deltax - pi, j * deltax - pi, deltax)) / a(j * deltax - pi));
		double aright = 1.0 / (theta(0, j * deltax - pi, (j + 1) * deltax - pi, deltax) / a(j * deltax - pi) + (1.0 - theta(0, j * deltax - pi, (j + 1) * deltax - pi, deltax)) / a((j + 1) * deltax - pi));
		u[i % 2][j] = lamuda * aleft*u[(i + 1) % 2][j - 1] + (1 - lamuda * aright - lamuda * aleft)*u[(i + 1) % 2][j] + lamuda * aright*u[(i + 1) % 2][j + 1];
	}
	u[i % 2][J] = 0;//�ұ߽紦��
}
//���ģ�����
double error_max(int J, double lamuda, double T, vector<vector<double>> u)
{
	double deltax = 2.0 * pi / J;              //�ռ���ɢ�̶ȣ�
	double deltat = lamuda * deltax * deltax;//���ռ���ɢ�̶ȣ�
	int N = floor(T / deltat) + 1;           //��Ҫǰ���Ĳ�����
	double max = 0;                 //��¼��
	for (int j = 0; j <= J; j++)
		if (fabs(u_exact(j*deltax - pi, T) - u[N % 2][j]) > max)
			max = fabs(u_exact(j*deltax - pi, T) - u[N % 2][j]);
	return max;
}
//���L2ģ��
double error_L2(int J, double lamuda, double T, vector<vector<double>> u)
{
	double deltax = 2.0 * pi / J;              //�ռ���ɢ�̶ȣ�
	double deltat = lamuda * deltax * deltax;//���ռ���ɢ�̶ȣ�
	int N = floor(T / deltat) + 1;           //��Ҫǰ���Ĳ�����
	double max = 0;                 //��¼��
	for (int j = 0; j <= J; j++)
		max += pow(u_exact(j*deltax - pi, T) - u[N % 2][j], 2);
	return sqrt(max*deltax);
}
//���������
void output_matrix(vector<vector<double>> matrix)
{
	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix[0].size(); j++)
			cout << matrix[i][j] << "\t";
		cout << endl;
	}
}
//�����ֵ�����ָ���ı��ļ���
void output(int i, double deltax, double lamuda, double T, vector<vector<double>> u)
{
	double deltat = lamuda * deltax;//���ռ���ɢ�̶ȣ�
	int N = T / deltat;             //��Ҫǰ���Ĳ�����
	ofstream outfile;//�����ά��ֵ��
	string temp = "E:\\study_materials\\PDEns\\06\\ns";
	string path = temp + to_string(i) + ".txt";
	outfile.open(path, ios::out);
	for (int j = 0; j < size(u[0]); j++)
		outfile << j * deltax << " " << u[N % 2][j] << endl;
}

int main()
{
	vector<vector<double>> error(5, vector<double>(5, 0));
	vector<vector<double>> error2(5, vector<double>(5, 0));
	for (int i = 1; i <= 5; i++)
	{
		int J = 10 * pow(2, i) + 1;                           //�ռ���ɢ�̶ȣ�
		double deltax = 2 * pi / J;                           
		double T = 1.0;                                       //����ʱ��T����ֵ�⣻
		vector<vector<double>> u(3, vector<double>(J + 1, 0));//��¼���������
		vector<vector<double>> u2(3, vector<double>(J + 1, 0));//��¼���������
		FTCS_mean(J, miu, T, u);
		FTCS_harmonic(J, miu, T, u2);
		error[i - 1][0] = J, error2[i - 1][0] = J;
		error[i - 1][1] = error_max(J, miu, T, u), error2[i - 1][1] = error_max(J, miu, T, u2);
		error[i - 1][3] = error_L2(J, miu, T, u), error2[i - 1][3] = error_L2(J, miu, T, u2);
	}
	for (int i = 1; i < 5; i++)
		error[i][2] = log(error[i - 1][1] / error[i][1]) / log((10 * pow(2, i + 1) + 1) / (10 * pow(2, i) + 1)),
		error[i][4] = log(error[i - 1][3] / error[i][3]) / log((10 * pow(2, i + 1) + 1) / (10 * pow(2, i) + 1)),
		error2[i][2] = log(error2[i - 1][1] / error2[i][1]) / log((10 * pow(2, i + 1) + 1) / (10 * pow(2, i) + 1)),
		error2[i][4] = log(error2[i - 1][3] / error2[i][3]) / log((10 * pow(2, i + 1) + 1) / (10 * pow(2, i) + 1));
	cout << "����ƽ����" << endl;
	cout << "N\tģ������\t������\tģL2���\t������" << endl;
	output_matrix(error);
	cout << endl << "����ƽ����" << endl;
	cout << "N\tģ������\t������\tģL2���\t������" << endl;
	output_matrix(error2);
}