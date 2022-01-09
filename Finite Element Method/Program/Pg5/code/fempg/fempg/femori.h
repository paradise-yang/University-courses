#pragma once
#include "improveCG.h"
#include <fstream>
#include <string>

using namespace std;

//真实解函数u(x)；
double u(double x)
{
	return (x - 1) * sin(x);
}
//函数f(x)；
double f(double x)
{
	return -2 * cos(x) + (x - 1) * sin(x);
}
//基于一次函数基生成系数矩阵K;
vector<vector<double>> make_cof_k(int N)
{
	//小区间个数有N个，系数矩阵cof_k大小为(N-1)*(N-1)；
	vector<vector<double>> cof_k(N - 1, vector<double>(N - 1, 0));
	for (int i = 1; i < N; i++)
		cof_k[i - 1][i - 1] = 2.0*N;
	for (int i = 0; i < N - 2; i++)
		cof_k[i][i + 1] = -1.0*N, cof_k[i + 1][i] = cof_k[i][i + 1];
	return cof_k;
	/*
	cof_k[0][0] = a[1][1] / he;
	for (int e = 2; e < N; e++)
	{
		cof_k[e - 2][e - 2] += (a[0][0] / he);
		cof_k[e - 1][e - 1] += (a[1][1] / he);
		cof_k[e - 2][e - 1] += (a[0][1] / he);
		cof_k[e - 1][e - 2] += (a[1][0] / he);
	}
	cof_k[N - 2][N - 2] += a[0][0] / he;
	*/
}
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
//基于一次函数的复化三点Gauss积分公式；
double gauss_integrate(int N, int i, double left, double right)
{
	int n = 4;
	double delta_x = (right - left) / n;
	double temp = 0;
	for (int j = 0; j < n; j++)
	{
		double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
		double x1 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 - sqrt(15) / 10 * delta_x;
		double x2 = (left + j * delta_x + left + (j + 1) * delta_x) / 2;
		double x3 = (left + j * delta_x + left + (j + 1) * delta_x) / 2 + sqrt(15) / 10 * delta_x;
		temp += a1 * f(x1) * basis_1(N, i, x1) + a2 * f(x2) * basis_1(N, i, x2) + a3 * f(x3) * basis_1(N, i, x3);
	}
	return temp;
}
//生成右侧向量；
vector<double> make_cof_f(int N)
{
	vector<double> cof_f(N - 1, 0);
	double he = 1.0 / N;
	for (int i = 1; i < N; i++)
		cof_f[i - 1] = gauss_integrate(N, i, (i - 1.0) * he, (i + 1.0) * he);
	return cof_f;
}
//基于1次基函数的求解函数uh(x)；
double uh_basis(int N, vector<double> cof_u, double x)
{
	double temp = 0;
	double he = 1.0 / N;
	for (int i = 1; i < N; i++)
		temp += cof_u[i - 1] * basis_1(N, i, x);
	return temp;
}

//离散计算无穷范数误差；
double error_max(int N, vector<double> cof_u)
{
	int n = 100;
	double delta_x = 1.0 / n;
	vector<double> error(n + 1, 0);
	for (int i = 0; i <= n; i++)
		error[i] = fabs(u(i * delta_x) - uh_basis(N, cof_u, i * delta_x));
	double max = error[0];
	for (int i = 1; i <= n; i++)
		if (error[i] > max)
			max = error[i];
	return max;
}
//基于复化Gauss积分公式离散计算L1范数误差；
double error_L1(int N, vector<double> cof_u)
{
	int n = 512;
	double delta_x = 1.0 / n;
	double temp = 0;
	for (int j = 0; j < n; j++)
	{
		double a1 = 5.0 / 18 * delta_x, a2 = 4.0 / 9 * delta_x, a3 = 5.0 / 18 * delta_x;
		double x1 = (j * delta_x + (j + 1) * delta_x) / 2 - sqrt(15) / 10 * delta_x;
		double x2 = (j * delta_x + (j + 1) * delta_x) / 2;
		double x3 = (j * delta_x + (j + 1) * delta_x) / 2 + sqrt(15) / 10 * delta_x;
		temp += a1 * fabs(u(x1) - uh_basis(N, cof_u, x1)) + a2 * fabs(u(x2) - uh_basis(N, cof_u, x2)) + a3 * fabs(u(x3) - uh_basis(N, cof_u, x3));
	}
	return temp;
}

//基函数计算求解，输出无穷范数误差；
void CG_error_out(vector<vector<double>>& error, int n)
{
	ofstream outfile;
	outfile.open("E:\\study_materials\\Finite Element Method\\Pg5\\code\\cg.txt", ios::out);
	for (int i = 0; i < n; i++)
	{
		int N = pow(2, i) * 10;                      //区间离散程度；
		error[i][0] = N;
		vector<vector<double>> cof_k = make_cof_k(N);//系数矩阵K；
		vector<double> cof_f = make_cof_f(N);        //右侧函数f与基内积向量；
		vector<double> cof_u(N - 1, 0);              //离散解的系数向量；
		int iteration_times = ITERATION_MAX;         //共轭梯度迭代法迭代次数；
		cof_u=conjugate_gradient(cof_k, cof_f,cof_u,PRECISION, iteration_times, error[i][1]);//共轭梯度法求解；
		//if (i == 0) matrix_output(cof_k), vector_output(cof_f), vector_output(cof_u);
		error[i][2] = error_max(N, cof_u);           //计算无穷范数作为误差；
		error[i][4] = error_L1(N, cof_u);            //计算L1范数作为误差；
		outfile << pow(2, i) * 10 << " " << error[i][1] << endl;//输入到指定文件；
		if (i > 0)
			error[i][3] = log(error[i - 1][2] / error[i][2]) / log(2),
			error[i][5] = log(error[i - 1][4] / error[i][4]) / log(2);
	}
	outfile.close();
}

//////以下为多重网格法求解
vector<double> multigrid_throw_up(vector<double> u)
{
	int dim_down = u.size(), dim_up = 2 * u.size() + 1;
	vector<double> up(dim_up, 0);
	for (int i = 0; i < dim_down; i++) up[2 * i + 1] = u[i];
	up[0] = u[0] / 2, up[dim_up - 1] = u[dim_down - 1] / 2;
	for (int i = 0; i < dim_down - 1; i++) up[2 * i + 2] = (u[i] + u[i + 1]) / 2;
	return up;
}
vector<double> multigrid_throw_down(vector<double> u)
{
	int dim_up = u.size(), dim_down = (u.size() - 1) / 2;
	vector<double> down(dim_down, 0);
	for (int i = 0; i < dim_down; i++)
		//down[i] = u[2 * i + 1],
		down[i] = (u[2 * i] + 2 * u[2 * i + 1] + u[2 * i + 2]) / 2;
	return down;
}

//生成第k层所需矩阵；
vector<vector<double>> Ak(int N)
{
	vector<vector<double>> cof_k(N - 1, vector<double>(N - 1, 0));
	for (int i = 1; i < N; i++) cof_k[i - 1][i - 1] = 2.0*N;
	for (int i = 0; i < N - 2; i++)	cof_k[i][i + 1] = -1.0*N, cof_k[i + 1][i] = cof_k[i][i + 1];
	return cof_k;
}

//迭代；
vector<double> iteration(vector<vector<double>> cof_Ak, vector<double> z, vector<double> g, double gamma_k)
{
	int dim = z.size();
	vector<double> temp = vector_subtract_vector(g, matrix_times_vector(cof_Ak, z));
	for (int i = 0; i < dim; i++)
		z[i] = z[i] + temp[i] / gamma_k;
	return z;
}
void iteration_two(int k, vector<double> &z, vector<double> g, int m1, double gammak)
{
	for (int t = 1; t <= m1; t++)
	{
		int i = 0;
		double tempi = z[0], temp = 0;
		z[0] = z[0] + g[0] / gammak - (2 * z[0] - z[1])*k / gammak;
		for (i = 1; i < z.size() - 1; i++)
		{
			temp = z[i];
			z[i] = z[i] + g[i] / gammak - (-tempi + 2 * z[i] - z[i + 1])*k / gammak;
			tempi = temp;
		}
		z[i] = z[i] + g[i] / gammak - (-tempi + 2 * z[i])*k / gammak;
	}
	/*
	int dim = z.size();
	vector<double> temp = vector_subtract_vector(g, matrix_times_vector(cof_Ak, z));
	for (int i = 0; i < dim; i++)
		z[i] = z[i] + temp[i] / gamma_k;
	return z;
	*/
}

//k-level multigrid kernel
vector<double> multigrid_kernel(int k_level, vector<double> z, vector<double> g, int m1, int m2, int p)
{
	int size0 = z.size(), n = z.size() + 1;
	vector<vector<double>> cof_Ak = Ak(n);
	if (k_level == 0)
	{
		int iteration_times = ITERATION_MAX;
		double time = 0;
		return conjugate_gradient(cof_Ak, g, z, PRECISION, iteration_times, time);
	}
	//presmooth:更新z,1——m1；
	double gamma_k = 4.0*n;
	for (int i = 1; i <= m1; i++) z = iteration(cof_Ak, z, g, gamma_k);
	//correction:更新z,m1——m1+1；
	vector<double> g_bar = multigrid_throw_down(vector_subtract_vector(g, matrix_times_vector(cof_Ak, z)));
	vector<double> q(g_bar);
	for (int i = 1; i <= p; i++) q = multigrid_kernel(k_level - 1, q, g_bar, m1, m2, p);
	z = vector_add_vector(z, multigrid_throw_up(q));
	//postsmooth:更新z,m1+1——m1+m2+1；
	for (int i = 1; i <= m2; i++) z = iteration(cof_Ak, z, g, gamma_k);
	return z;
}

//k-level multigrid kernel
vector<double> multigrid_kernel_two(int k_level, vector<double> z, vector<double> g, int m1, int m2, int p)
{
	int size0 = z.size(), n = z.size() + 1;
	vector<vector<double>> cof_Ak = Ak(n);
	if (k_level == 0)
	{
		int iteration_times = ITERATION_MAX;
		double time = 0;
		return conjugate_gradient(cof_Ak, g, z, PRECISION, iteration_times, time);
	}
	//presmooth:更新z,1——m1；
	double gamma_k = 4.0*n;
	iteration_two(n, z, g, m1, gamma_k);
	//for (int i = 1; i <= m1; i++) z = iteration(cof_Ak, z, g, gamma_k);
	//correction:更新z,m1——m1+1；
	int templow = (size0 - 1) / 2;
	vector<double> g_bar(templow, 0);
	int i = 0;
	g_bar[0] = (g[0] + 2 * g[1] + g[2] - (2 * z[1] - z[3])*n) / 2;
	for (i = 1; i < templow - 1; i++)
		g_bar[i] = (g[2 * i] + 2 * g[2 * i + 1] + g[2 * i + 2]
			- (-z[2 * i - 1] + 2 * z[2 * i + 1] - z[2 * i + 3])*n) / 2;
	g_bar[i] = (g[2 * i] + 2 * g[2 * i + 1] + g[2 * i + 2] - (-z[2 * i - 1] + 2 * z[2 * i + 1])*n) / 2;
	vector<double> q(templow, 0);
	//vector<double> g_bar = multigrid_throw_down(vector_subtract_vector(g, matrix_times_vector(cof_Ak, z)));
	//vector<double> q(g_bar);
	for (int i = 1; i <= p; i++) q = multigrid_kernel_two(k_level - 1, q, g_bar, m1, m2, p);
	for (int j = 0; j < templow; j++) z[2 * j + 1] += q[j];
	z[0] += q[0] / 2, z[z.size() - 1] += q[templow - 1] / 2;
	for (int j = 0; j < templow - 1; j++) z[2 * j + 2] += (q[j] + q[j + 1]) / 2;
	//z = vector_add_vector(z, multigrid_throw_up(q));
	//postsmooth:更新z,m1+1——m1+m2+1；
	iteration_two(n, z, g, m2, gamma_k);
	//for (int i = 1; i <= m2; i++) z = iteration(cof_Ak, z, g, gamma_k);
	return z;
}

//full Multigrid  Algorithm
vector<double> full_multigrid(int k_level, int r, int m1, int m2, int p, double& time)
{
	/*
	if (k_level == 0)
	{
		int size0 = u0.size(), n = u0.size() + 1;
		vector<vector<double>> cof_k = make_cof_k(n);
		vector<double> cof_f = make_cof_f(n);
		int iteration_times = ITERATION_MAX;
		return conjugate_gradient(cof_k, cof_f, u0, PRECISION, iteration_times, time);
	}
	int size0 = 2 * u0.size() + 1, n = 2 * u0.size() + 2;
	vector<double> cof_f = make_cof_f(n);
	vector<double> u = multigrid_throw_up(u0);//上一步进行上投影；
	auto start = std::chrono::steady_clock::now();
	for (int i = 1; i <= r; i++) u = multigrid_kernel(k_level, u, cof_f, m1, m2, p);//r步迭代；
	auto end = std::chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	time = elapsed_seconds.count();
	return u;
	*/
	int n = 10;
	vector<vector<double>> cof_k = make_cof_k(n);
	vector<double> cof_f = make_cof_f(n);
	vector<double> cof_u(n - 1, 0);
	auto start = std::chrono::steady_clock::now();
	int iteration_times = ITERATION_MAX;
	vector<double> u = conjugate_gradient(cof_k, cof_f, cof_u, PRECISION, iteration_times, time);
	if (k_level == 0)
		return u;
	for (int level = 1; level <= k_level; level++)
	{
		u = multigrid_throw_up(u);
		vector<double> cof_g = make_cof_f(n*pow(2, level));
		for (int i = 0; i < r; i++)
			u = multigrid_kernel(level, u, cof_g, m1, m2, p);
	}
	auto end = std::chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	time = elapsed_seconds.count();
	return u;
}

vector<double> full_multigrid_two(int k_level, int r, int m1, int m2, int p, double& time)
{
	int n = 10;
	vector<vector<double>> cof_k = make_cof_k(n);
	vector<double> cof_f = make_cof_f(n);
	vector<double> cof_u(n - 1, 0);
	auto start = std::chrono::steady_clock::now();
	int iteration_times = ITERATION_MAX;
	vector<double> u = conjugate_gradient(cof_k, cof_f, cof_u, PRECISION, iteration_times, time);
	if (k_level == 0)
		return u;
	for (int level = 1; level <= k_level; level++)
	{
		u = multigrid_throw_up(u);
		vector<double> cof_g = make_cof_f(n*pow(2, level));
		for (int i = 0; i < r; i++)
			u = multigrid_kernel_two(level, u, cof_g, m1, m2, p);
	}
	auto end = std::chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	time = elapsed_seconds.count();
	return u;
}

//基函数计算求解，输出无穷范数误差；
void MG_error_out(vector<vector<double>>& error, int n)
{
	ofstream outfile;
	outfile.open("E:\\study_materials\\Finite Element Method\\Pg5\\code\\mg.txt", ios::out);
	for (int i = 0; i < n; i++)
	{
		int N = pow(2, i) * 10;                             //区间离散程度；
		error[i][0] = N;
		int r = 100, m1 = 3, m2 = 3, p = 2;                   //相关参数；
		vector<double> u = full_multigrid(i, r, m1, m2, p, error[i][1]);
		//u = full_multigrid(i, u, r, m1, m2, p, error[i][1]);//多重网格求解；
		//if (i == 0) vector_output(u);
		error[i][2] = error_max(N, u);                      //计算无穷范数作为误差；
		error[i][4] = error_L1(N, u);                       //计算L1范数作为误差；
		outfile << pow(2, i) * 10 << " " << error[i][1] << endl;//输入到指定文件；
		if (i > 0)
			//error[i][1] += error[i - 1][1],
			error[i][3] = log(error[i - 1][2] / error[i][2]) / log(2),
			error[i][5] = log(error[i - 1][4] / error[i][4]) / log(2);
	}
	outfile.close();
}

void MG_error_out_two(vector<vector<double>>& error, int n)
{
	ofstream outfile;
	outfile.open("E:\\study_materials\\Finite Element Method\\Pg5\\code\\MG.txt", ios::out);
	for (int i = 0; i < n; i++)
	{
		int N = pow(2, i) * 10;                             //区间离散程度；
		error[i][0] = N;
		int r = 100, m1 = 3, m2 = 3, p = 2;                   //相关参数；
		vector<double> u = full_multigrid_two(i, r, m1, m2, p, error[i][1]);
		//u = full_multigrid(i, u, r, m1, m2, p, error[i][1]);//多重网格求解；
		//if (i == 0) vector_output(u);
		error[i][2] = error_max(N, u);                      //计算无穷范数作为误差；
		error[i][4] = error_L1(N, u);                       //计算L1范数作为误差；
		outfile << pow(2, i) * 10 << " " << error[i][1] << endl;//输入到指定文件；
		if (i > 0)
			//error[i][1] += error[i - 1][1],
			error[i][3] = log(error[i - 1][2] / error[i][2]) / log(2),
			error[i][5] = log(error[i - 1][4] / error[i][4]) / log(2);
	}
	outfile.close();
}

void CG_error_out_two(vector<vector<double>>& error, int n)
{
	ofstream outfile;
	outfile.open("E:\\study_materials\\Finite Element Method\\Pg5\\code\\cg2.txt", ios::out);
	for (int i = 0; i < n; i++)
	{
		int N = pow(2, i) * 10;                      //区间离散程度；
		error[i][0] = N;
		vector<vector<double>> cof_k = make_cof_k(N);//系数矩阵K；
		vector<double> cof_f = make_cof_f(N);        //右侧函数f与基内积向量；
		vector<double> cof_u(N - 1, 0);              //离散解的系数向量；
		int iteration_times = ITERATION_MAX;         //共轭梯度迭代法迭代次数；
		//cof_u = conjugate_gradient(cof_k, cof_f, cof_u, PRECISION, iteration_times, error[i][1]);//共轭梯度法求解；
		auto start_time = std::chrono::steady_clock::now();
		cof_u = Conjugate_Gradient(cof_k, cof_f, cof_u, PRECISION, iteration_times);
		auto end_time = std::chrono::steady_clock::now();
		chrono::duration<double> elapsed_seconds = end_time - start_time;
		error[i][1] = elapsed_seconds.count();
		//if (i == 0) matrix_output(cof_k), vector_output(cof_f), vector_output(cof_u);
		error[i][2] = error_max(N, cof_u);           //计算无穷范数作为误差；
		error[i][4] = error_L1(N, cof_u);            //计算L1范数作为误差；
		outfile << pow(2, i) * 10 << " " << error[i][1] << endl;//输入到指定文件；
		if (i > 0)
			error[i][3] = log(error[i - 1][2] / error[i][2]) / log(2),
			error[i][5] = log(error[i - 1][4] / error[i][4]) / log(2);
	}
	outfile.close();
}