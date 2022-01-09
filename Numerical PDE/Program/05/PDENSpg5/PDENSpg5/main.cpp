#include<iostream>
#include<iomanip>
#include<vector>
#include<cmath>
#include<fstream>
#include<string>

using namespace std;

#define pi 3.1415926535897932384626433832795028841971 //Բ����\piֵ��

//��ֵ������
double u_0(double x)
{
	if (x >= 0.4 && x <= 0.6)
		return 1.0;
	else
		return 0.0;
}
//��ȷ�⺯����
double u_exact(double x, double t)
{
	if ((x - t >= 0.4) && (x - t <= 0.6))
		return 1.0;
	else
		return 0.0;
}
//FTBS��ʽ��
void FTBS(double deltax, double lamuda, double T, vector<vector<double>> &u)
{
	int J = 7.0 / deltax;           //�ռ�������������
	double deltat = lamuda * deltax;//���ռ���ɢ�̶ȣ�
	int N = T / deltat;             //��Ҫǰ���Ĳ�����
	//��ֵ��
	for (int j = 0; j <= J; j++)
		u[0][j] = u_0(j*deltax);
	//FTBS��ʽ��
	for (int i = 1; i <= N; i++)
	{
		u[i % 2][0] = u[(i + 1) % 2][0] - lamuda * (u[(i + 1) % 2][0] - 0);//u[(i + 1) % 2][J - 1]);
		for (int j = 1; j <= J; j++)
			u[i % 2][j] = u[(i + 1) % 2][j] - lamuda * (u[(i + 1) % 2][j] - u[(i + 1) % 2][j - 1]);
	}
}
//���ģ�����
double error_max(double deltax, double lamuda, double T, vector<vector<double>> u)
{
	int J = 6.0 / deltax;           //�ռ�������������
	double deltat = lamuda * deltax;//���ռ���ɢ�̶ȣ�
	int N = T / deltat;             //��Ҫǰ���Ĳ�����
	double max = 0;                 //��¼��
	for (int j = 0; j <= J; j++)
		if (fabs(u_exact(j*deltax, T) - u[N % 2][j]) > max)
			max = fabs(sin(2 * pi * (j * deltax + T)) - u[N % 2][j]);
	return max;
}
//�����ֵ�����ָ���ı��ļ���
void output(int i, double deltax, double lamuda, double T, vector<vector<double>> u)
{
	double deltat = lamuda * deltax;//���ռ���ɢ�̶ȣ�
	int N = T / deltat;             //��Ҫǰ���Ĳ�����
	ofstream outfile;//�����ά��ֵ��
	string temp = "E:\\study_materials\\PDEns\\05\\ns";
	string path = temp + to_string(i) + ".txt";
	outfile.open(path, ios::out);
	for (int j = 0; j < size(u[0]); j++)
		outfile << j * deltax << " " << u[N % 2][j] << endl;
}

int main()
{
	double deltax = 0.05;                                  //�ռ���ɢ�̶ȣ�
	int J = 7.0 / deltax;
	double T = 1.0;                                        //����ʱ��T����ֵ�⣻
	double lamuda = 0.2;                                   //r=\lamuda=deltat/deltax��
	vector<vector<double>> u(3, vector<double>(J + 1, 0)); //��¼���������
	double max = 0;                                        //��¼��
	
	cout << endl << "deltax=0.05,lamuda=0.2,T=1.0��" << endl;
	FTBS(deltax, lamuda, T, u);
	max = error_max(deltax, lamuda, T, u);
	cout << "ģ������Ϊ��" << max << endl;
	//�����ֵ�⣻
	output(0, deltax, lamuda, T, u);
	cout << "��ֵ�����μ��ĵ�" << endl << endl;
	
	cout << "deltax=0.05,lamuda=0.2,T=2.0��" << endl;
	T = 2.0;
	FTBS(deltax, lamuda, T, u);
	max = error_max(deltax, lamuda, T, u);
	cout << "ģ������Ϊ��" << max << endl;
	//�����ֵ�⣻
	output(1, deltax, lamuda, T, u);
	cout << "��ֵ�����μ��ĵ�" << endl << endl;

	cout << "deltax=0.05,lamuda=0.2,T=5.0��" << endl;
	T = 5.0;
	FTBS(deltax, lamuda, T, u);
	max = error_max(deltax, lamuda, T, u);
	cout << "ģ������Ϊ��" << max << endl;
	//�����ֵ�⣻
	output(2, deltax, lamuda, T, u);
	cout << "��ֵ�����μ��ĵ�" << endl << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	cout << endl << "deltax=0.05,lamuda=0.8,T=1.0��" << endl;
	lamuda = 0.8;
	T = 1.0;
	FTBS(deltax, lamuda, T, u);
	max = error_max(deltax, lamuda, T, u);
	cout << "ģ������Ϊ��" << max << endl;
	//�����ֵ�⣻
	output(3, deltax, lamuda, T, u);
	cout << "��ֵ�����μ��ĵ�" << endl << endl;

	cout << "deltax=0.05,lamuda=0.8,T=2.0��" << endl;
	T = 2.0;
	FTBS(deltax, lamuda, T, u);
	max = error_max(deltax, lamuda, T, u);
	cout << "ģ������Ϊ��" << max << endl;
	//�����ֵ�⣻
	output(4, deltax, lamuda, T, u);
	cout << "��ֵ�����μ��ĵ�" << endl << endl;

	cout << "deltax=0.05,lamuda=0.8,T=5.0��" << endl;
	T = 5.0;
	FTBS(deltax, lamuda, T, u);
	max = error_max(deltax, lamuda, T, u);
	cout << "ģ������Ϊ��" << max << endl;
	//�����ֵ�⣻
	output(5, deltax, lamuda, T, u);
	cout << "��ֵ�����μ��ĵ�" << endl << endl;
}