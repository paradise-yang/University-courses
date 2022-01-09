#include "Iteration.h"

using namespace std;

typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

int main()
{
	//Ch4.1用Jacobi、Gauss-Seidel、SOR求解离散形式微分方程，epsilon=1.0 / 0.1 / 0.01 / 0.0001；
	
	double epsilon = 1.0;     //此处修改epsilon值；
	double w = 1.75;           //SOR迭代的松弛因子；
	double error = 1e-5;
	clock_t start_time, end_time;
	int dim1 = 100, times = 0;
	
	double a = 0.5, h = 1.0 / dim1;//保留4位有效数字，观察精确解，所以取
	double ah2 = a * h * h;
	Matrix matrix(dim1 - 1, Vector(dim1 - 1, 0.0)), matrix0;
	Vector x0(dim1 - 1, 1.0), value(dim1 - 1, ah2), exact(dim1), y(dim1 - 1, 0);
	value[dim1 - 2] = value[dim1 - 2] - epsilon - h;
	
	cout << endl << "Ch4.1用Jacobi、Gauss-Seidel、SOR求解离散形式微分方程，epsilon = " << epsilon << endl;
	F(dim1, epsilon, a, y);   //求解精确解；
	cout << "精确解为：" << endl;
	Vector_cout(y);

	cout << endl << "Jacobi迭代：";
	band_matrix(dim1 - 1, matrix, epsilon, h);
	matrix0 = matrix;
	exact = x0;
	start_time = clock();                //获取第一个时间；
	Jacobi_iteration(dim1 - 1, matrix, x0, value, error, times);
	end_time = clock();                  //获取第二个时间；
	cout << "迭代解为：" << endl;
	Vector_cout(x0);
	cout << "迭代次数为 = " << times << endl;
	cout << "运行时间为 =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	Vector_subtraction(x0, y, x0);
	cout << "与精确解误差||x - x*|| = " << vector_infinity_norm(x0) << endl;
	cout << "与精确解误差||x - x*||2 = " << vector_2_norm(x0) << endl;
	matrix.clear();
	x0.clear();

	cout << endl << "Gauss-Seidel迭代：";
	matrix = matrix0;
	x0 = exact;
	start_time = clock();                //获取第一个时间；
	Gauss_Seidel_iteration(dim1 - 1, matrix, x0, value, error, times);
	end_time = clock();                  //获取第二个时间；
	cout << "迭代解为：" << endl;
	Vector_cout(x0);
	cout << "迭代次数为 = " << times << endl;
	cout << "运行时间为 =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	Vector_subtraction(x0, y, x0);
	cout << "与精确解误差||x - x*|| = " << vector_infinity_norm(x0) << endl;
	cout << "与精确解误差||x - x*||2 = " << vector_2_norm(x0) << endl;
	matrix.clear();
	x0.clear();

	cout << endl << "SOR迭代：" << endl << "松弛因子为 w = " << w << endl;
	matrix = matrix0;
	x0 = exact;
	start_time = clock();                //获取第一个时间；
	SOR_iteration(dim1 - 1, w, matrix, x0, value, error, times);
	end_time = clock();                  //获取第二个时间；
	cout << "迭代解为：" << endl;
	Vector_cout(x0);
	cout << "迭代次数为 = " << times << endl;
	cout << "运行时间为 =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	Vector_subtraction(x0, y, x0);
	cout << "与精确解误差||x - x*|| = " << vector_infinity_norm(x0) << endl;
	cout << "与精确解误差||x - x*||2 = " << vector_2_norm(x0) << endl;
	matrix.clear();
	x0.clear();
	
	//Ch4.2 Jacobi、Gauss-Seidel、SOR迭代求解偏微分方程，N=20 / 40 /80；
	error = 1e-7;
	cout << endl << "Ch4.2 Guass-Seidel迭代求解偏微分方程，N=20 / 40 /80：" << endl;
	
	int N = 20;
	cout << "N = 20 ：" << endl;
	
	start_time = clock();                //获取第一个时间；
	Jacobi_PDE(N, times, error);
	end_time = clock();                  //获取第二个时间；
	cout << "Jacobi迭代次数 = " << times;
	cout << "           运行时间为 =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	
	start_time = clock();                //获取第一个时间；
	Gauss_Seidel_PDE(N, times, error);
	end_time = clock();                  //获取第二个时间；
	cout << "Gauss-Seidel迭代次数 = " << times;
	cout << "     运行时间为 =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;

	start_time = clock();                //获取第一个时间；
	SOR_PDE(N, times, error, w);
	end_time = clock();                  //获取第二个时间；
	cout << "SOR迭代次数 = " << times;
	cout << "               运行时间为 =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl << endl;
	
	N = 40;
	cout << "N = 40 ：" << endl;

	start_time = clock();                //获取第一个时间；
	Jacobi_PDE(N, times, error);
	end_time = clock();                  //获取第二个时间；
	cout << "Jacobi迭代次数 = " << times;
	cout << "          运行时间为 =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;

	start_time = clock();                //获取第一个时间；
	Gauss_Seidel_PDE(N, times, error);
	end_time = clock();                  //获取第二个时间；
	cout << "Gauss-Seidel迭代次数 = " << times;
	cout << "     运行时间为 =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;

	start_time = clock();                //获取第一个时间；
	SOR_PDE(N, times, error, w);
	end_time = clock();                  //获取第二个时间；
	cout << "SOR迭代次数 = " << times;
	cout << "              运行时间为 =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl << endl;
	
	N = 80;
	cout << "N = 80 ：" << endl;

	start_time = clock();                //获取第一个时间；
	Jacobi_PDE(N, times, error);
	end_time = clock();                  //获取第二个时间；
	cout << "Jacobi迭代次数 = " << times;
	cout << "          运行时间为 =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;

	start_time = clock();                //获取第一个时间；
	Gauss_Seidel_PDE(N, times, error);
	end_time = clock();                  //获取第二个时间；
	cout << "Gauss-Seidel迭代次数 = " << times;
	cout << "    运行时间为 =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;

	start_time = clock();                //获取第一个时间；
	SOR_PDE(N, times, error, w);
	end_time = clock();                  //获取第二个时间；
	cout << "SOR迭代次数 = " << times;
	cout << "              运行时间为 =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl << endl;
	
}