#include "Conjugate_gradient_method.h"

using namespace std;

typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

int main() 
{
	int dim = 20;
	int times = 0;
	double w = 1.725;            //SOR迭代的松弛因子；
	double error = 1e-7;         //SOR迭代的精度；
	double epsilon = 1e-10;      //实用共轭梯度迭代法的epsiloon；
	double iteration_error = 0;  //共轭梯度迭代法的输出参考误差；
	clock_t start_time, end_time;
	Matrix matrix((dim - 1)* (dim - 1), Vector((dim - 1)* (dim - 1), 0));
	Vector b((dim - 1) * (dim - 1), 0), x((dim - 1) * (dim - 1), 0);
	

	//Ch5.1 共轭梯度法求解PDE，并与SOR方法比较，寻找最佳松弛因子；
	cout << "Ch5.1 共轭梯度法求解PDE，并与SOR方法比较，寻找最佳松弛因子。" << endl;
	Ch5_1(dim, matrix, b);
	Conjugate_Gradient((dim - 1) * (dim - 1), matrix, b, epsilon, x, times, iteration_error);
	cout << endl << "共轭梯度法迭代次数为 = " << times << endl;
	cout << "参考误差为 = " << iteration_error << endl;
	cout << "迭代解为：" << endl;
	Vector_cout(x);
	
	cout << endl;
	Vector x0((dim - 1) * (dim - 1), 1.0);
	start_time = clock();                //获取第一个时间；
	SOR_iteration((dim - 1) * (dim - 1), w, matrix, x0, b, error, times);
	//Vector_cout(x0);
	end_time = clock();                  //获取第二个时间；
	cout << "SOR迭代法松弛因子为w = " << w << endl;
	cout << "迭代次数为 = " << times << endl;
	cout << "运行时间为 =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	


	//Ch5.2 用Hilbert矩阵测试共轭梯度法；
	cout << endl << endl <<endl << "Ch5.2 用Hilbert矩阵测试共轭梯度法。" << endl;
	
	dim = 100;     //Hilbert矩阵维数；
	Matrix hilbert(dim,Vector(dim,0));
	Vector value(dim, 0), y(dim, 0);
	Hilbert_matrix(dim, hilbert);
	for (int i = 0; i < dim; i++)                //初始化值；
	{
		for (int j = 0; j < dim; j++)
			value[i] += hilbert[i][j];
		value[i] = value[i] / 3.0;
	}
	start_time = clock();                //获取第一个时间；
	Conjugate_Gradient(dim, hilbert, value, epsilon, y, times, iteration_error);
	end_time = clock();                  //获取第二个时间；
	cout << endl << dim << "阶Hilbert方阵：" << "迭代次数：" << times << endl;
	cout << "运行时间为 =  " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	cout << "参考误差为 = " << iteration_error << endl;
	cout << "迭代解为：" << endl;
	Vector_cout(y);
	


	//Ch5.3分别用Jacobi、Gauss-Seidel、共轭梯度法求解方程组；
	cout << endl << endl << endl << "Ch5.3分别用Jacobi、Gauss-Seidel、共轭梯度法求解方程组。" << endl;
	dim = 5;
	Matrix A = { {10.0,1.0,2.0,3.0,4.0},{1.0,9.0,-1.0,2.0,-3.0},{2.0,-1.0,7.0,3.0,-5.0},{3.0,2.0,3.0,12.0,-1.0},{4.0,-3.0,-5.0,-1.0,15.0} };
	Vector z = { 12.0,-27.0,14.0,-17.0,12.0 };
	Vector z0(dim, 1.0);
	
	cout << endl << "Jacobi迭代：" << endl;
	Jacobi_iteration(dim, A, z0, z, error, times);
	cout << "迭代次数为 = " << times << endl;
	cout << "迭代解为：";
	Vector_cout(z0);
	
	A = { {10.0,1.0,2.0,3.0,4.0},{1.0,9.0,-1.0,2.0,-3.0},{2.0,-1.0,7.0,3.0,-5.0},{3.0,2.0,3.0,12.0,-1.0},{4.0,-3.0,-5.0,-1.0,15.0} };
	z0 = { 1.0,1.0,1.0,1.0,1.0 };
	cout << endl << "Gauss-Seidel迭代：" << endl;
	Gauss_Seidel_iteration(dim, A, z0, z, error, times);
	cout << "迭代次数为 = " << times << endl;
	cout << "迭代解为：";
	Vector_cout(z0);

	A = { {10.0,1.0,2.0,3.0,4.0},{1.0,9.0,-1.0,2.0,-3.0},{2.0,-1.0,7.0,3.0,-5.0},{3.0,2.0,3.0,12.0,-1.0},{4.0,-3.0,-5.0,-1.0,15.0} };
	z0 = { 1.0,1.0,1.0,1.0,1.0 };
	cout << endl << "共轭梯度法：" << endl;
	Conjugate_Gradient(dim, A, z, epsilon, z0, times, iteration_error);
	cout << "迭代次数为 = " << times << endl;
	cout << "参考误差为 = " << iteration_error << endl;
	cout << "迭代解为：";
	Vector_cout(z0);
	
}