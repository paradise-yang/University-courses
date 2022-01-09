#pragma once
#include"femori.h"

using namespace std;

//一次有限元空间///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Shinshkin下构造一次函数基；
double basis_shinshkin_1(int N, int i, double left, double right, double x)
{
	double he = (right - left) / N;
	if ((x > left + (i - 1.0) * he) && (x <= left + i * he))
		return (x - left - (i - 1.0) * he) / he;
	else if ((x > left + i * he) && (x < left + (i + 1.0) * he))
		return (left + (i + 1.0) * he - x) / he;
	else
		return 0.0;
}
//Shinshkin下构造一次函数基的导数；
double basis_shinshkin_diff_1(int N, int i, double left, double right, double x)
{
	double he = (right - left) / N;
	if ((x > left + (i - 1.0) * he) && (x <= left + i * he))
		return 1.0 / he;
	else if ((x > left + i * he) && (x < left + (i + 1.0) * he))
		return -1.0 / he;
	else
		return 0.0;
}
//基于1次基函数的求解函数uh(x)；
double uh_basis_shinshkin_1(int N, double tao, vector<double> cof_u, double x)
{
	double he_1 = tao / N;
	double he_2 = (-tao + 1.0) / N;
	double temp = 0;
	for (int i = 1; i < N; i++)
		//cout<<endl<< cof_u_1[i - 1] * basis_shinshkin_1(N, i, 0, tao, x) + cof_u_2[i - 1] * basis_shinshkin_1(N, i, tao, 1 - tao, x),
		temp += (cof_u[i - 1] * basis_shinshkin_1(N, i, 0, tao, x) + cof_u[i - 1 + N] * basis_shinshkin_1(N, i, tao, 1, x));

	int i = N;
	if ((x > (i - 1.0) * he_1) && (x <= tao))
		return temp + cof_u[i - 1] * (x - (i - 1.0) * he_1) / he_1;
	else if ((x > tao) && (x < tao + he_2))
		return temp + cof_u[i - 1] * (tao + he_2 - x) / he_2;
	else
		return temp;
}

//Shinshkin下构造系数矩阵
double gauss_integrate_shinshkin_k(int N, double basisleft, double basisright, int i, int k, double left, double right)
{
	int n = 2;
	double delta_x = (right - left) / n;
	double temp = 0;
	for (int j = 0; j < n; j++)
	{
		double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
		double x1 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 - sqrt(15) / 10 * delta_x;
		double x2 = (left + j * delta_x + left + (j + 1) * delta_x) / 2;
		double x3 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 + sqrt(15) / 10 * delta_x;
		temp += a1 * (basis_shinshkin_diff_1(N, i, basisleft, basisright, x1) * basis_shinshkin_diff_1(N, k, basisleft, basisright, x1) * epsilon7 + basis_shinshkin_1(N, i, basisleft, basisright, x1) * basis_shinshkin_diff_1(N, k, basisleft, basisright, x1))
			+ a2 * (basis_shinshkin_diff_1(N, i, basisleft, basisright, x2) * basis_shinshkin_diff_1(N, k, basisleft, basisright, x2) * epsilon7 + basis_shinshkin_1(N, i, basisleft, basisright, x2) * basis_shinshkin_diff_1(N, k, basisleft, basisright, x2))
			+ a3 * (basis_shinshkin_diff_1(N, i, basisleft, basisright, x3) * basis_shinshkin_diff_1(N, k, basisleft, basisright, x3) * epsilon7 + basis_shinshkin_1(N, i, basisleft, basisright, x3) * basis_shinshkin_diff_1(N, k, basisleft, basisright, x3));
	}
	return temp;
}
//Shinshkin下基于一次函数基生成系数矩阵K;
vector<vector<double>> make_cof_shinshkin_k(int N)
{
	double tao = -2 * epsilon7 * log(N) + 1.0;
	double he_1 = tao / N;//离散区间长度；
	double he_2 = (1.0 - tao) / N;
	vector<vector<double>> cof_k(2 * N - 1, vector<double>(2 * N - 1, 0));
	for (int e = 1; e < N - 1; e++)
	{
		cof_k[e - 1][e - 1] = 2 * (epsilon7 / he_1 + 1.0);//(epsilon + he_1) * 2 / he_1;
		cof_k[e - 1][e] = -(epsilon7 / he_1 + 1.0) + 0.5;//-(epsilon + he_1) / he_1 + 0.5;
		cof_k[e][e - 1] = -(epsilon7 / he_1 + 1.0) - 0.5;//-(epsilon + he_1) / he_1 - 0.5;

		cof_k[e - 1 + N][e - 1 + N] = 2 * (epsilon7 / he_2 + 1.0);//(epsilon + he_2) * 2 / he_2;
		cof_k[e - 1 + N][e - 1 + N + 1] = -(epsilon7 / he_2 + 1.0) + 0.5;//-(epsilon + he_2) / he_2 + 0.5;
		cof_k[e - 1 + N + 1][e - 1 + N] = -(epsilon7 / he_2 + 1.0) - 0.5;//-(epsilon + he_2) / he_2 - 0.5;
	}
	int e = N - 1;
	cof_k[e - 1][e - 1] = 2 * (epsilon7 / he_1 + 1.0);//(epsilon + he_1) * 2 / he_1;
	cof_k[e - 1][e] = -(epsilon7 / he_1 + 1.0) + 0.5;//-(epsilon + he_1) / he_1 + 0.5;
	cof_k[e][e - 1] = -(epsilon7 / he_1 + 1.0) - 0.5;//-(epsilon + he_1) / he_1 - 0.5;

	cof_k[e - 1 + N][e - 1 + N] = 2 * (epsilon7 / he_2 + 1.0); //(epsilon + he_2) * 2 / he_2;

	e = N;//第n个区间长度为he_1
	cof_k[e - 1][e - 1] = epsilon7 / he_1 + epsilon7 / he_2 + 2.0;//(epsilon + he_1) / he_1 + (epsilon + he_2) / he_2;
	cof_k[e - 1][e] = -(epsilon7 / he_2 + 1.0) + 0.5; //-(epsilon + he_2) / he_2 + 0.5;
	cof_k[e][e - 1] = -(epsilon7 / he_2 + 1.0) - 0.5; //-(epsilon + he_2) / he_2 - 0.5;

	return cof_k;
}

//Shinshkin下基于一次函数的复化三点Gauss积分公式；
double gauss_integrate_shinshkin_1(int N, double basisleft, double basisright, int i, double left, double right)
{
	int n = 2;
	double delta_x = (right - left) / n;
	double temp = 0;
	for (int j = 0; j < n; j++)
	{
		double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
		double x1 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 - sqrt(15) / 10 * delta_x;
		double x2 = (left + j * delta_x + left + (j + 1) * delta_x) / 2;
		double x3 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 + sqrt(15) / 10 * delta_x;
		temp += a1 * f(x1) * basis_shinshkin_1(N, i, basisleft, basisright, x1) + a2 * f(x2) * basis_shinshkin_1(N, i, basisleft, basisright, x2) + a3 * f(x3) * basis_shinshkin_1(N, i, basisleft, basisright, x3);
	}
	return temp;
}
//生成一次空间的右侧向量；
void make_cof_shinshkin_f(int N, double tao, vector<double>& cof_f)
{
	double he_1 = tao / N;
	double he_2 = (1.0 - tao) / N;
	for (int i = 1; i < N; i++)
	{
		cof_f[i - 1] = (i - 1.0) * he_1 * he_1;
		cof_f[i - 1 + N] = ((i - 1.0) * he_2 + tao) * he_2;// + he_2 * tao;
		//cof_f_1[i - 1] = gauss_integrate_shinshkin_1(N, 0, tao, i, (i - 1) * he_1, (i + 1) * he_1);
		//cof_f_2[i - 1] = gauss_integrate_shinshkin_1(N, tao, 1 - tao, i, tao+(i - 1) * he_2, tao+(i + 1) * he_2);
	}
	int i = N;
	cof_f[i - 1] = -he_1 * he_1*2.0 / 3 - he_2 * he_2 / 3 + 3.0 * tao*he_1 / 2 - tao * he_2 / 2;//he_1 / 2 + (-2 * tao + 1.0) / 6 / N / N;
}

double error_max_shinshkin_1(int N, double tao, vector<double> cof_u)
{
	ofstream outfile;
	string path = "E:\\study_materials\\Finite Element Method\\Program\\Pg6\\output\\shinshkin17" + to_string(N) + ".txt";
	outfile.open(path, ios::out);
	int n = 100000;
	double delta_x = 1.0 / n;
	vector<double> error(n + 1, 0);
	vector<double> u(n - 1, 0);
	for (int i = 0; i <= n; i++)
		outfile << i * delta_x << " " << uh_basis_shinshkin_1(N, tao, cof_u, i * delta_x) << endl,
		error[i] = fabs(u_2(i * delta_x) - uh_basis_shinshkin_1(N, tao, cof_u, i * delta_x));
	outfile.close();
	double max = error[0];
	for (int i = 1; i <= n; i++)
		if (error[i] > max)
			max = error[i];
	return max;
}
void error_L1_L2_shinshkin_1(int N, double tao, vector<double> cof_u, double& errorL1, double& errorL2)
{
	int n = 1024;
	//double delta_x = 1.0 / n;
	double deltax1 = tao / n, deltax2 = (1.0 - tao) / n;
	double temp1 = 0, temp11 = 0, temp2 = 0, temp22 = 0;
	double df1 = 0, df2 = 0, df3 = 0, dff1 = 0, dff2 = 0, dff3 = 0;
	double a1 = 5.0 / 18 * deltax1, a2 = 4.0 / 9 * deltax1, a3 = 5.0 / 18 * deltax1;
	double aa1 = 5.0 / 18 * deltax2, aa2 = 4.0 / 9 * deltax2, aa3 = 5.0 / 18 * deltax2;
	for (int j = 0; j < n; j++)
	{
		double x1 = (j * deltax1 + (j + 1) * deltax1) / 2 - sqrt(15) / 10 * deltax1;
		double x2 = (j * deltax1 + (j + 1) * deltax1) / 2;
		double x3 = (j * deltax1 + (j + 1) * deltax1) / 2 + sqrt(15) / 10 * deltax1;
		df1 = fabs(u_2(x1) - uh_basis_shinshkin_1(N, tao, cof_u, x1));
		df2 = fabs(u_2(x2) - uh_basis_shinshkin_1(N, tao, cof_u, x2));
		df3 = fabs(u_2(x3) - uh_basis_shinshkin_1(N, tao, cof_u, x3));
		temp1 += a1 * df1 + a2 * df2 + a3 * df3;
		temp2 += a1 * df1*df1 + a2 * df2*df2 + a3 * df3*df3;
		double xx1 = (2 * tao + j * deltax1 + (j + 1) * deltax1) / 2 - sqrt(15) / 10 * deltax1;
		double xx2 = (2 * tao + j * deltax1 + (j + 1) * deltax1) / 2;
		double xx3 = (2 * tao + j * deltax1 + (j + 1) * deltax1) / 2 + sqrt(15) / 10 * deltax1;
		dff1 = fabs(u_2(x1) - uh_basis_shinshkin_1(N, tao, cof_u, xx1));
		dff2 = fabs(u_2(x2) - uh_basis_shinshkin_1(N, tao, cof_u, xx2));
		dff3 = fabs(u_2(x3) - uh_basis_shinshkin_1(N, tao, cof_u, xx3));
		temp11 += aa1 * dff1 + aa2 * dff2 + aa3 * dff3;
		temp22 += aa1 * dff1* dff1 + aa2 * dff2* dff2 + aa3 * dff3* dff3;
	}
	errorL1 = temp1 + temp11;
	errorL2 = sqrt(temp2 + temp22);
}
//基函数计算求解，输出无穷范数误差；
void shinshkin_error_out(int n, vector<vector<double>>& error, int flag)
{
	ofstream outfile;
	string path = "E:\\study_materials\\Finite Element Method\\Program\\Pg6\\shinshkin1" + to_string(flag) + ".txt";
	outfile.open(path, ios::out);
	for (int i = 0; i < n; i++)
	{
		//区间离散程度；
		int N = pow(2, i) * 10;
		error[i][0] = N;
		double tao = 1.0 - 2 * epsilon7 * log(N);
		vector<vector<double>> cof_k_shinshkin = make_cof_shinshkin_k(N);
		vector<double> cof_f_shinshkin(2 * N - 1, 0);
		make_cof_shinshkin_f(N, tao, cof_f_shinshkin);
		vector<double> cof_u_shinshkin(2 * N - 1, 0);
		//if (i == 0) matrix_output(cof_k_shinshkin), cout << endl, vector_output(cof_f_shinshkin), cout << endl, vector_output(cof_u_shinshkin), cout << endl;
		//int times = 1000;             //共轭梯度迭代法迭代次数；
		//cof_u_shinshkin = conjugate_gradient_quick(cof_k_shinshkin, cof_f_shinshkin, cof_u_shinshkin, 1e-20, times, error[i][1]);
		//if (i == 0)	matrix_output(cof_k_shinshkin), cout << endl, vector_output(cof_f_shinshkin);//check运算输出；
		//auto start_time = std::chrono::steady_clock::now();
		//cof_u_shinshkin = gauss_line(cof_k_shinshkin, cof_f_shinshkin);
		gauss(cof_k_shinshkin, cof_f_shinshkin, cof_u_shinshkin);
		//auto end_time = chrono::steady_clock::now();
		//chrono::duration<double> elapsed_seconds = end_time - start_time;
		//error[i][1] = elapsed_seconds.count();
		//if (i == 0)	cout << endl, vector_output(cof_u_shinshkin);//check运算输出；
		error[i][2] = error_max_shinshkin_1(N, tao, cof_u_shinshkin);
		error_L1_L2_shinshkin_1(N, tao, cof_u_shinshkin, error[i][4], error[i][6]);
		//error[i][4] = error_L1_shinshkin_1(N, tao, cof_u_shinshkin);
		if (i > 0)
			error[i][3] = log(error[i - 1][2] / error[i][2]) / log(2),
			error[i][5] = log(error[i - 1][4] / error[i][4]) / log(2),
			error[i][7] = log(error[i - 1][6] / error[i][6]) / log(2);
		outfile << N << " " << error[i][2] << " " << error[i][3] << " " << error[i][4] << " " << error[i][5] << " " << error[i][6] << " " << error[i][7] << " " << endl;//输入到指定文件；
	}
	outfile.close();
}



//一次有限元空间///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<vector<double>> make_cof_shinshkin_k_2(int N)
{
	double tao = -2 * epsilon7 * log(N) + 1.0;
	double he1 = tao / N;//离散区间长度；
	double he2 = (1.0 - tao) / N;
	vector<vector<double>> cof_k(4 * N - 1, vector<double>(4 * N - 1, 0));
	for (int i = 0; i < N - 1; i++)
	{
		cof_k[i][i] = (epsilon7 / he1 + 1.0)*14.0 / 3;
		cof_k[i][i + N - 1] = -(epsilon7 / he1 + 1.0) * 8.0 / 3 - 2.0 / 3, cof_k[i][i + N] = -(epsilon7 / he1 + 1.0) * 8.0 / 3 + 2.0 / 3;
		cof_k[i + N - 1][i] = -(epsilon7 / he1 + 1.0) * 8.0 / 3 + 2.0 / 3, cof_k[i + N][i] = -(epsilon7 / he1 + 1.0) * 8.0 / 3 - 2.0 / 3;

		cof_k[i + 2 * N][i + 2 * N] = (epsilon7 / he2 + 1.0) * 14.0 / 3;
		cof_k[i + 2 * N][i + 3 * N - 1] = -(epsilon7 / he2 + 1.0) * 8.0 / 3 - 2.0 / 3, cof_k[i + 2 * N][i + 3 * N] = -(epsilon7 / he2 + 1.0) * 8.0 / 3 + 2.0 / 3;
		cof_k[i + 3 * N - 1][i + 2 * N] = -(epsilon7 / he2 + 1.0) * 8.0 / 3 + 2.0 / 3, cof_k[i + 3 * N][i + 2 * N] = -(epsilon7 / he2 + 1.0) * 8.0 / 3 - 2.0 / 3;
	}
	for (int i = 0; i < N - 2; i++)
	{
		cof_k[i][i + 1] = (epsilon7 / he1 + 1.0) / 3 - 1.0 / 6, cof_k[i + 1][i] = (epsilon7 / he1 + 1.0) / 3 + 1.0 / 6;
		cof_k[i + 2 * N][i + 2 * N + 1] = (epsilon7 / he2 + 1.0) / 3 - 1.0 / 6, cof_k[i + 2 * N + 1][i + 2 * N] = (epsilon7 / he2 + 1.0) / 3 + 1.0 / 6;
	}
	for (int i = 0; i < N; i++)
	{
		cof_k[i + N - 1][i + N - 1] = (epsilon7 / he1 + 1.0) * 16.0 / 3;
		cof_k[i + 3 * N - 1][i + 3 * N - 1] = (epsilon7 / he2 + 1.0) * 16.0 / 3;
	}
	int e = 2 * N - 1;
	cof_k[e][e] = (epsilon7 / he1 + 1.0)*7.0 / 3 + (epsilon7 / he2 + 1.0)*7.0 / 3;
	cof_k[e][e + 1] = (epsilon7 / he2 + 1.0) / 3 - 1.0 / 6;
	cof_k[e + 1][e] = (epsilon7 / he2 + 1.0) / 3 + 1.0 / 6;
	cof_k[e][e - 1] = -(epsilon7 / he1 + 1.0) * 8.0 / 3 - 2.0 / 3;
	cof_k[e - 1][e] = -(epsilon7 / he1 + 1.0) * 8.0 / 3 + 2.0 / 3;
	cof_k[e][N - 2] = (epsilon7 / he1 + 1.0) / 3 + 1.0 / 6;
	cof_k[N - 2][e] = (epsilon7 / he1 + 1.0) / 3 - 1.0 / 6;
	cof_k[e][3 * N - 1] = -(epsilon7 / he2 + 1.0) * 8.0 / 3 + 2.0 / 3;
	cof_k[3 * N - 1][e] = -(epsilon7 / he2 + 1.0) * 8.0 / 3 - 2.0 / 3;

	return cof_k;
}
vector<double> make_cof_shinshkin_f_2(int N,double tao)
{
	vector<double> cof_f(4 * N - 1, 0);
	double he1 = tao / N;//离散区间长度；
	double he2 = (1.0 - tao) / N;
	for (int i = 1; i < N; i++)
		cof_f[i - 1] = (i - 1)*1.0 / 3 * he1*he1,
		cof_f[2 * N + i - 1] = (i - 1)*1.0 / 3 * he2*he2;
	for (int i = 0; i < N; i++)
		cof_f[i + N - 1] = (2 * i - 1)*1.0 / 3 * he1*he1,
		cof_f[3 * N + i - 1] = (2 * i - 1)*1.0 / 3 * he2*he2;

	int i = 2 * N - 1;
	cof_f[i] = -he1 * he1 / 6 - he2 * he2 / 6 + 7.0*he1*tao / 6 - he2 * tao*5.0 / 6;

	return cof_f;
}

double basis_shinshkin_2(int N, int i, double tao, double x)
{
	double he1 = tao / N;//离散区间长度；
	double he2 = (1.0 - tao) / N;
	if (i < N)
	{
		if (x > (i - 1) * he1 && x <= i * he1)
			return (2 * x - (2 * i - 1) * he1) * (x - (i - 1) * he1) / he1 / he1;
		else if (x > i * he1 && x < (i + 1) * he1)
			return (2 * x - (2 * i + 1) * he1) * (x - (i + 1) * he1) / he1 / he1;
		else
			return 0;
	}
	else if(i<2*N)
	{
		if (x > (i - N) * he1 && x < (i - N + 1) * he1)
			return 4 * (x - (i - N) * he1) * ((i - N + 1) * he1 - x) / he1 / he1;
		else
			return 0;
	}
	else if (i < 3 * N)
	{
		int e = i - 2 * N;
		if (x > (tao + (e - 1)*he2) && x <= (tao + e * he2))
			return (2 * x - 2 * tao - (2 * e - 1) * he2) * (x - tao - (e - 1) *he2) / he2 / he2;
		else if (x > (tao + e * he2) && x < (tao + (e + 1)*he2))
			return (2 * x - 2 * tao - (2 * i + 1) * he2) * (x - tao - (i + 1) * he2) / he2 / he2;
		else
			return 0.0;
	}
	else
	{
		int e = i - 2 * N;
		if (x > (tao + (e - N)*he2) && x < (tao + (e - N + 1)*he2))
			return 4 * (x - tao - (e - N) * he2) * ((e - N + 1) * he2 + tao - x) / he2 / he2;
		else
			return 0.0;
	}
}
double uh_basis_shinshkin_2(int N, double tao, vector<double> cof_u, double x)
{
	int dim = 2 * N - 1;
	double he_1 = tao / N;
	double he_2 = (-tao + 1.0) / N;
	double temp = 0;
	if (x < tao)
	{
		for (int i = 0; i <= dim; i++)
			temp += cof_u[i] * basis_shinshkin_2(N, i + 1, tao, x);
		return temp;
	}
	else
	{
		int size = cof_u.size();
		for (int i = dim; i < size; i++)
			temp += cof_u[i] * basis_shinshkin_2(N, i + 1, tao, x);
		return temp;
	}
}

double error_max_shinshkin_2(int N, double tao, vector<double> cof_u)
{
	ofstream outfile;
	string path = "E:\\study_materials\\Finite Element Method\\Program\\Pg6\\output\\shinshkin27" + to_string(N) + ".txt";
	outfile.open(path, ios::out);
	int n = 100000;
	double delta_x = 1.0 / n;
	vector<double> error(n + 1, 0);
	vector<double> u(n - 1, 0);
	for (int i = 0; i <= n; i++)
		outfile << i * delta_x << " " << uh_basis_shinshkin_2(N, tao, cof_u, i * delta_x) << endl,
		error[i] = fabs(u_2(i * delta_x) - uh_basis_shinshkin_2(N, tao, cof_u, i * delta_x));
	outfile.close();
	double max = error[0];
	for (int i = 1; i <= n; i++)
		if (error[i] > max)
			max = error[i];
	return max;
}
void error_L1_L2_shinshkin_2(int N, double tao, vector<double> cof_u, double& errorL1, double& errorL2)
{
	int n = 1024;
	//double delta_x = 1.0 / n;
	double deltax1 = tao / n, deltax2 = (1.0 - tao) / n;
	double temp1 = 0, temp11 = 0, temp2 = 0, temp22 = 0;
	double df1 = 0, df2 = 0, df3 = 0, dff1 = 0, dff2 = 0, dff3 = 0;
	double a1 = 5.0 / 18 * deltax1, a2 = 4.0 / 9 * deltax1, a3 = 5.0 / 18 * deltax1;
	double aa1 = 5.0 / 18 * deltax2, aa2 = 4.0 / 9 * deltax2, aa3 = 5.0 / 18 * deltax2;
	for (int j = 0; j < n; j++)
	{
		double x1 = (j * deltax1 + (j + 1) * deltax1) / 2 - sqrt(15) / 10 * deltax1;
		double x2 = (j * deltax1 + (j + 1) * deltax1) / 2;
		double x3 = (j * deltax1 + (j + 1) * deltax1) / 2 + sqrt(15) / 10 * deltax1;
		df1 = fabs(u_2(x1) - uh_basis_shinshkin_2(N, tao, cof_u, x1));
		df2 = fabs(u_2(x2) - uh_basis_shinshkin_2(N, tao, cof_u, x2));
		df3 = fabs(u_2(x3) - uh_basis_shinshkin_2(N, tao, cof_u, x3));
		temp1 += a1 * df1 + a2 * df2 + a3 * df3;
		temp2 += a1 * df1*df1 + a2 * df2*df2 + a3 * df3*df3;
		double xx1 = (2 * tao + j * deltax1 + (j + 1) * deltax1) / 2 - sqrt(15) / 10 * deltax1;
		double xx2 = (2 * tao + j * deltax1 + (j + 1) * deltax1) / 2;
		double xx3 = (2 * tao + j * deltax1 + (j + 1) * deltax1) / 2 + sqrt(15) / 10 * deltax1;
		dff1 = fabs(u_2(x1) - uh_basis_shinshkin_2(N, tao, cof_u, xx1));
		dff2 = fabs(u_2(x2) - uh_basis_shinshkin_2(N, tao, cof_u, xx2));
		dff3 = fabs(u_2(x3) - uh_basis_shinshkin_2(N, tao, cof_u, xx3));
		temp11 += aa1 * dff1 + aa2 * dff2 + aa3 * dff3;
		temp22 += aa1 * dff1* dff1 + aa2 * dff2* dff2 + aa3 * dff3* dff3;
	}
	errorL1 = temp1 + temp11;
	errorL2 = sqrt(temp2 + temp22);
}
void shinshkin_error_out_2(int n, vector<vector<double>>& error, int flag)
{
	ofstream outfile;
	string path = "E:\\study_materials\\Finite Element Method\\Program\\Pg6\\shinshkin2" + to_string(flag) + ".txt";
	outfile.open(path, ios::out);
	for (int i = 0; i < n; i++)
	{
		//区间离散程度；
		int N = pow(2, i) * 10;
		error[i][0] = N;
		double tao = 1.0 - 2 * epsilon7 * log(N);
		vector<vector<double>> cof_k_shinshkin = make_cof_shinshkin_k_2(N);
		vector<double> cof_f_shinshkin = make_cof_shinshkin_f_2(N, tao);
		vector<double> cof_u_shinshkin(4 * N - 1, 0);
		//auto start_time = std::chrono::steady_clock::now();
		cof_u_shinshkin = gauss_line(cof_k_shinshkin, cof_f_shinshkin);
		//cof_u_shinshkin = cholesky(cof_k_shinshkin, cof_f_shinshkin);
		//auto end_time = chrono::steady_clock::now();
		//chrono::duration<double> elapsed_seconds = end_time - start_time;
		//error[i][1] = elapsed_seconds.count();
		//if (i == 0)	cout << endl, vector_output(cof_u_shinshkin);//check运算输出；
		error[i][2] = error_max_shinshkin_2(N, tao, cof_u_shinshkin);
		error_L1_L2_shinshkin_2(N, tao, cof_u_shinshkin, error[i][4], error[i][6]);//计算L1范数误差与L2范数误差；
		//error[i][4] = error_L1_shinshkin_2(N, tao, cof_u_shinshkin);
		if (i > 0)
			error[i][3] = log(error[i - 1][2] / error[i][2]) / log(2),
			error[i][5] = log(error[i - 1][4] / error[i][4]) / log(2),
			error[i][7] = log(error[i - 1][6] / error[i][6]) / log(2);
		outfile << N << " " << error[i][2] << " " << error[i][3] << " " << error[i][4] << " " << error[i][5] << " " << error[i][6] << " " << error[i][7] << " " << endl;//输入到指定文件；
	}
	outfile.close();
}
