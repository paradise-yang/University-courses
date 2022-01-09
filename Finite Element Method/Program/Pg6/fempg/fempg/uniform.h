#pragma once
#include"femori.h"

using namespace std;

//均匀划分，一次有限元空间//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//第i个一次基函数；
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
//基于1次基函数的求解函数uh(x)；
double uh_basis_1(int N, vector<double> cof_u, double x)
{
	double temp = 0;
	double he = 1.0 / N;
	for (int i = 1; i < N; i++)
		temp += cof_u[i - 1] * basis_1(N, i, x);
	return temp;
}

//基于一次函数基生成系数矩阵K;
vector<vector<double>> make_cof_k(int N, double epsilon)
{
	//离散区间有N个，系数矩阵cof_k大小为(N-1)*(N-1)；//double he = 1.0 / N;
	vector<vector<double>> cof_k(N - 1, vector<double>(N - 1, 0));
	for (int i = 0; i < N - 1; i++)
		cof_k[i][i] = 2.0*(epsilon*N + 1.0);
	for (int i = 0; i < N - 2; i++)
		cof_k[i][i + 1] = 0.5 - (epsilon*N + 1.0),
		cof_k[i + 1][i] = -0.5 - (epsilon*N + 1.0);
	return cof_k;
	/*
	for (int e = 1; e < N - 1; e++)
	{
		cof_k[e - 1][e - 1] = 2 * (epsilon*N + 1);//gauss_integrate_k(N, e, e, (e - 1) * he, (e + 1) * he);

		cof_k[e - 1][e] = 0.5 - (epsilon*N + 1);//gauss_integrate_k(N, e, e + 1, e * he, (e + 1) * he);
		cof_k[e][e - 1] = 0.5 + (epsilon*N + 1);//gauss_integrate_k(N, e + 1, e, e * he, (e + 1) * he);
	}
	int e = N - 1;
	cof_k[e - 1][e - 1] = gauss_integrate_k(N, e, e, (e - 1) * he, (e + 1) * he);

	cof_k[0][0] = (7.0 / 3 + he * he / 5 + (1 - cos(he)) / he) / he;
	for (int e = 2; e < N; e++)
	{
		double a = (e - 1.0) * he;
		cof_k[e - 2][e - 2] += ((7.0 / 3 + a * a / 3 + a * he / 6 + he * he / 30 + (cos(a) - cos(a + he)) / he) / he);
		cof_k[e - 1][e - 1] += ((7.0 / 3 + a * a / 3 + a * he / 2 + he * he / 5 + (cos(a) - cos(a + he)) / he) / he);
		cof_k[e - 2][e - 1] += ((-11.0 / 6 + a * a / 6 + a * he / 6 + he * he / 20 + (-cos(a) + cos(a + he)) / he) / he);
		cof_k[e - 1][e - 2] += ((-11.0 / 6 + a * a / 6 + a * he / 6 + he * he / 20 + (-cos(a) + cos(a + he)) / he) / he);
	}
	double a = (N - 1.0) * he;
	cof_k[N - 2][N - 2] += ((7.0 / 3 + a * a / 3 + a * he / 6 + he * he / 30 + (cos(a) - cos(a + he)) / he) / he);
	*/
}
//生成一次空间的右侧向量；
vector<double> make_cof_f(int N)
{
	double he = 1.0 / N;
	vector<double> cof_f(N - 1, 0);
	for (int i = 1; i < N; i++)
		cof_f[i - 1] = (i - 1)*1.0 / N / N;// *he*he; // N / N;//gauss_integrate_1(N, i, (i - 1) * he, (i + 1) * he);
	return cof_f;
}

double error_max(int N, vector<double> cof_u, double epsilon)
{
	ofstream outfile;
	string path = "E:\\study_materials\\Finite Element Method\\Program\\Pg6\\output\\uniform17" + to_string(N) + ".txt";
	outfile.open(path, ios::out);
	int n = 50000;
	double delta_x = 1.0 / n;
	vector<double> error(n + 1, 0);
	for (int i = 0; i <= n; i++)
		outfile << i * delta_x << " " << uh_basis_1(N, cof_u, i * delta_x) << endl,
		error[i] = fabs(u(epsilon, i * delta_x) - uh_basis_1(N, cof_u, i * delta_x));
	outfile.close();
	double max = error[0];
	for (int i = 1; i <= n; i++)
		if (error[i] > max)
			max = error[i];
	return max;
}
void error_L1_L2(int N, vector<double> cof_u, double epsilon, double& errorL1, double& errorL2)
{
	int n = 256;
	double delta_x = 1.0 / n;
	double temp1 = 0, temp2 = 0, df1 = 0, df2 = 0, df3 = 0, df4 = 0;
	for (int j = 0; j < n; j++)
	{
		double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
		double x1 = (j * delta_x + (j + 1) * delta_x) / 2 - sqrt(15) / 10 * delta_x;
		double x2 = (j * delta_x + (j + 1) * delta_x) / 2;
		double x3 = (j * delta_x + (j + 1) * delta_x) / 2 + sqrt(15) / 10 * delta_x;
		df1 = fabs(u(epsilon, x1) - uh_basis_1(N, cof_u, x1));
		df2 = fabs(u(epsilon, x2) - uh_basis_1(N, cof_u, x2));
		df3 = fabs(u(epsilon, x3) - uh_basis_1(N, cof_u, x3));
		temp1 += a1 * df1 + a2 * df2 + a3 * df3;
		temp2 += a1 * df1*df1 + a2 * df2*df2 + a3 * df3*df3;
		//temp += a1 * fabs(u(epsilon, x1) - uh_basis_1(N, cof_u, x1)) + a2 * fabs(u(epsilon, x2) - uh_basis_1(N, cof_u, x2)) + a3 * fabs(u(epsilon, x3) - uh_basis_1(N, cof_u, x3));
	}
	errorL1 = temp1;
	errorL2 = sqrt(temp2);
}
void CG_error_out(vector<vector<double>>& error, int n, double epsilon, int flag)
{
	ofstream outfile;
	string path = "E:\\study_materials\\Finite Element Method\\Program\\Pg6\\uniform1" + to_string(flag) + ".txt";
	outfile.open(path, ios::out);
	for (int i = 0; i < n; i++)
	{
		int N = pow(2, i) * 10;                      //区间离散程度；
		error[i][0] = N;
		vector<vector<double>> cof_k = make_cof_k(N, epsilon);//系数矩阵K；
		vector<double> cof_f = make_cof_f(N);        //右侧函数f与基内积向量；
		vector<double> cof_u(N - 1, 1);              //离散解的系数向量；
		//int iteration_times = ITERATION_MAX;         //共轭梯度迭代法迭代次数；
		//cof_u = conjugate_gradient_quick(cof_k, cof_f, cof_u, PRECISION, iteration_times, error[i][1]);//共轭梯度法求解；
		//gauss(cof_k, cof_f, cof_u);
		cof_u = gauss_line(cof_k, cof_f);
		//if (i == 0) matrix_output(cof_k), vector_output(cof_f), vector_output(cof_u);
		error[i][2] = error_max(N, cof_u,epsilon);               //计算无穷范数作为误差；
		error_L1_L2(N, cof_u, epsilon, error[i][4], error[i][6]);//计算L1范数误差与L2范数误差；
		//error[i][4] = error_L1(N, cof_u,epsilon);            //计算L1范数作为误差；
		if (i > 0)
			error[i][3] = log(error[i - 1][2] / error[i][2]) / log(2),
			error[i][5] = log(error[i - 1][4] / error[i][4]) / log(2),
			error[i][7] = log(error[i - 1][6] / error[i][6]) / log(2);
		outfile << N << " " << error[i][2] << " " << error[i][3] << " " << error[i][4] << " " << error[i][5] << " " << error[i][6] << " " << error[i][7] << " " << endl;//输入到指定文件；
	}
	outfile.close();
}



//均匀划分，二次有限元空间///////////////////////////////////////////////////////////////////////////////////////////////////////////

//第i个二次基函数；
double basis_2(int N, int i, double x)
{
	double he = 1.0 / N;
	if (i < N)
	{
		if (x > (i - 1) * he && x <= i * he)
			return (2 * x - (2 * i - 1) * he) * (x - (i - 1) * he) / he / he;
		else if (x > i * he && x < (i + 1) * he)
			return (2 * x - (2 * i + 1) * he) * (x - (i + 1) * he) / he / he;
		else
			return 0;
	}
	else
	{
		if (x > (i - N) * he && x < (i - N + 1) * he)
			return 4 * (x - (i - N) * he) * ((i - N + 1) * he - x) / he / he;
		else
			return 0;
	}
}
//基于2次基函数的求解函数uh(x)；
double uh_basis_2(int N, vector<double> cof_u, double x)
{
	double temp = 0.0;
	double he = 1.0 / N;
	for (int i = 1; i < 2 * N; i++)
		temp += cof_u[i - 1] * basis_2(N, i, x);
	return temp;
}

//基于二次函数基生成系数矩阵K；
vector<vector<double>> make_cof_k_2(int N, double epsilon)
{
	//小区间个数有N个，系数矩阵cof_k大小为(2N-1)*(2N-1)；
	double he = 1.0 / N;//离散区间长度；
	vector<vector<double>> cof_k(2 * N - 1, vector<double>(2 * N - 1, 0));
	for (int i = 0; i < N - 1; i++)	
		cof_k[i][i] = 14.0 / 3 * (epsilon*N + 1);
	for (int i = 0; i < N - 2; i++)
		cof_k[i][i + 1] = 1.0 / 3 * (epsilon*N + 1) - 1.0 / 6,
		cof_k[i + 1][i] = 1.0 / 3 * (epsilon*N + 1) + 1.0 / 6;
	for (int i = 0; i < N - 1; i++)
		cof_k[i][i + N - 1] = -8.0 / 3 * (epsilon*N + 1) - 2.0 / 3, cof_k[i][i + N] = -8.0 / 3 * (epsilon*N + 1) + 2.0 / 3,
		cof_k[i + N - 1][i] = -8.0 / 3 * (epsilon*N + 1) + 2.0 / 3, cof_k[i + N][i] = -8.0 / 3 * (epsilon*N + 1) - 2.0 / 3;
	for (int i = 0; i < N; i++) 
		cof_k[i + N - 1][i + N - 1] = 16.0 / 3 * (epsilon*N + 1);
	/*
	for (int i = 1; i < N - 1; i++)
	{
		cof_k[i - 1][i - 1] = 14.0 / 3 * N;//gauss_integrate_k_2(N, i, i, (i - 1) * he, (i + 1) * he);

		cof_k[i - 1][i] = gauss_integrate_k_2(N, i, i + 1, i * he, (i + 1) * he);
		cof_k[i][i - 1] = cof_k[i - 1][i];

		cof_k[i - 1][i - 1 + N - 1] = gauss_integrate_k_2(N, i, i + N - 1, (i - 1) * he, i * he);
		cof_k[i - 1 + N - 1][i - 1] = cof_k[i - 1][i - 1 + N - 1];

		cof_k[i - 1][i - 1 + N] = gauss_integrate_k_2(N, i, i + N, i * he, (i + 1) * he);
		cof_k[i - 1 + N][i - 1] = cof_k[i - 1][i - 1 + N];

		cof_k[i - 1 + N - 1][i - 1 + N - 1] = gauss_integrate_k_2(N, i + N - 1, i + N - 1, (i - 1) * he, i * he);
	}
	int i = N - 1;
	cof_k[N - 2][N - 2] = gauss_integrate_k_2(N, i, i, (i - 1) * he, (i + 1) * he);

	cof_k[N - 1 - 1][N - 1 - 1 + N - 1] = gauss_integrate_k_2(N, i, i + N - 1, (i - 1) * he, i * he);
	cof_k[N - 1 - 1 + N - 1][N - 1 - 1] = cof_k[N - 1 - 1][N - 1 - 1 + N - 1];

	cof_k[N - 1 - 1][N - 1 - 1 + N] = gauss_integrate_k_2(N, i, i + N, i * he, (i + 1) * he);
	cof_k[N - 1 - 1 + N][N - 1 - 1] = cof_k[N - 1 - 1][N - 1 - 1 + N];

	cof_k[2 * N - 3][2 * N - 3] = gauss_integrate_k_2(N, i + N - 1, i + N - 1, (i - 1) * he, i * he);

	i = N;
	cof_k[2 * N - 2][2 * N - 2] = gauss_integrate_k_2(N, i + N - 1, i + N - 1, (i - 1) * he, i * he);
	*/
	return cof_k;
}
//生成右侧向量；
vector<double> make_cof_f_2(int N)
{
	double he = 1.0 / N;
	vector<double> cof_f(2 * N - 1, 0);
	for (int i = 1; i < N; i++)
		cof_f[i - 1] = (i - 1)*1.0 / 3 / N / N;//gauss_integrate_2(N, i, (i - 1) * he, (i + 1) * he);
	for (int i = 0; i < N; i++)
		cof_f[i + N - 1] = (2 * i - 1)*1.0 / 3 / N / N;
	//for (int i = N; i < 2 * N; i++)
	//	cof_f[i - 1] = gauss_integrate_2(N, i, (i - N) * he, (i - N + 1) * he);
	return cof_f;
}

double error_max_2(int N, vector<double> cof_u, double epsilon)
{
	ofstream outfile;
	string path = "E:\\study_materials\\Finite Element Method\\Program\\Pg6\\output\\uniform27" + to_string(N) + ".txt";
	outfile.open(path, ios::out);
	int n = 50000;
	double delta_x = 1.0 / n;
	vector<double> error(n + 1, 0);
	for (int i = 0; i <= n; i++)
		outfile << i * delta_x << " " << uh_basis_2(N, cof_u, i * delta_x) << endl,
		error[i] = fabs(u(epsilon, i * delta_x) - uh_basis_2(N, cof_u, i * delta_x));
	outfile.close();
	double max = error[0];
	for (int i = 1; i <= n; i++)
		if (error[i] > max)
			max = error[i];
	return max;
}
void error_L1_L2_2(int N, vector<double> cof_u, double epsilon, double& errorL1, double& errorL2)
{
	int n = 256;
	double delta_x = 1.0 / n;
	double temp1 = 0, temp2 = 0, df1 = 0, df2 = 0, df3 = 0, df4 = 0;
	for (int j = 0; j < n; j++)
	{
		double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
		double x1 = (j * delta_x + (j + 1) * delta_x) / 2 - sqrt(15) / 10 * delta_x;
		double x2 = (j * delta_x + (j + 1) * delta_x) / 2;
		double x3 = (j * delta_x + (j + 1) * delta_x) / 2 + sqrt(15) / 10 * delta_x;
		df1 = fabs(u(epsilon, x1) - uh_basis_2(N, cof_u, x1));
		df2 = fabs(u(epsilon, x2) - uh_basis_2(N, cof_u, x2));
		df3 = fabs(u(epsilon, x3) - uh_basis_2(N, cof_u, x3));
		temp1 += a1 * df1 + a2 * df2 + a3 * df3;
		temp2 += a1 * df1*df1 + a2 * df2*df2 + a3 * df3*df3;
		//temp += a1 * fabs(u(epsilon, x1) - uh_basis_2(N, cof_u, x1)) + a2 * fabs(u(epsilon, x2) - uh_basis_2(N, cof_u, x2)) + a3 * fabs(u(epsilon, x3) - uh_basis_2(N, cof_u, x3));
	}
	errorL1 = temp1;
	errorL2 = sqrt(temp2);
}
void CG_error_out_2(vector<vector<double>>& error, int n, double epsilon, int flag)
{
	ofstream outfile;
	string path = "E:\\study_materials\\Finite Element Method\\Program\\Pg6\\uniform2"+ to_string(flag) +".txt";
	outfile.open(path, ios::out);
	for (int i = 0; i < n; i++)
	{
		int N = pow(2, i) * 10;                      //区间离散程度；
		error[i][0] = N;
		vector<vector<double>> cof_k = make_cof_k_2(N,epsilon);//系数矩阵K；
		vector<double> cof_f = make_cof_f_2(N);        //右侧函数f与基内积向量；
		vector<double> cof_u(2 * N - 1, 0.0);            //离散解的系数向量；
		//int iteration_times = ITERATION_MAX;           //共轭梯度迭代法迭代次数；
		//cof_u = conjugate_gradient(cof_k, cof_f, cof_u, 1e-20, iteration_times, error[i][1]);//共轭梯度法求解；
		//gauss(cof_k, cof_f, cof_u);
		//cof_u = cholesky(cof_k, cof_f);
		cof_u = gauss_line(cof_k, cof_f);
		//if (i == 0)
		//	matrix_output(cof_k), cout << endl,
		//	vector_output(cof_f), cout << endl,
		//	vector_output(cof_u);
		error[i][2] = error_max_2(N, cof_u,epsilon);             //计算无穷范数作为误差；
		error_L1_L2_2(N, cof_u, epsilon, error[i][4], error[i][6]);//计算L1范数误差与L2范数误差；
		//error[i][4] = error_L1_2(N, cof_u,epsilon);            //计算L1范数作为误差；
		if (i > 0)
			error[i][3] = log(error[i - 1][2] / error[i][2]) / log(2),
			error[i][5] = log(error[i - 1][4] / error[i][4]) / log(2),
			error[i][7] = log(error[i - 1][6] / error[i][6]) / log(2);
		outfile << N << " " << error[i][2] << " " << error[i][3] << " " << error[i][4] << " " << error[i][5] << " " << error[i][6] << " " << error[i][7] << endl;//输入到指定文件；
	}
	outfile.close();
}