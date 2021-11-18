#include "Symmetric_eigenvalue.h"

using namespace std;

typedef vector<vector<double>> Matrix;
typedef vector<double> Vector;

int main()
{
	//Ch7.1用Jacobi过关法求50-100阶对称三角阵的全部特征值；
	cout << "Ch7.1 用Jacobi过关法求50-100阶对称三角阵的全部特征值。" << endl << endl;
	cout << "请问，想输出某个阶数方阵特征值，还是50-100阶所有均输出？ 前者输入1，后者输入0：";
	int temp = 0;
	double u = 1e-7;//精度，也即关值deta；
	cin >> temp;
	if (temp == 1)
	{
		int dim = 0;
		cout << endl << "请输入想查看的矩阵阶数：";
		cin >> dim;
		Matrix matrix(dim, Vector(dim, 0));
		for (int i = 0; i < dim; i++) matrix[i][i] = 4;
		for (int i = 1; i < dim; i++) matrix[i][i - 1] = 1, matrix[i - 1][i] = 1;
		Pass_Jacobi_method(matrix, u);
		Matrix_cout(matrix, dim, dim);
	}
	else
	{
		for (int dim = 50; dim <= 100; dim++)
		{
			Matrix matrix(dim, Vector(dim, 0));
			for (int i = 0; i < dim; i++) matrix[i][i] = 4;
			for (int i = 1; i < dim; i++) matrix[i][i - 1] = 1, matrix[i - 1][i] = 1;
			Pass_Jacobi_method(matrix, u);
		}
	}



	//Ch7.2用二分法求实对称三角阵指定特征值；
	cout << endl << endl << endl;
	cout << "Ch7.2 用二分法求实对称三角阵指定特征值。" << endl;
	int dim = 100;
	Matrix matrix(dim, Vector(dim, 0));
	for (int i = 0; i < dim; i++) matrix[i][i] = 2;
	for (int i = 1; i < dim; i++) matrix[i][i - 1] = -1, matrix[i - 1][i] = -1;
	cout << "矩阵第一个特征值（最小特征值）为 " << bisection_method(matrix, 1) << endl;
	cout << "矩阵第100个特征值（最大特征值）为 " << bisection_method(matrix, dim) << endl;
}