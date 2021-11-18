#include "Norm_Accuracy.h"

using namespace std;

typedef  vector<vector<double>> Matrix;

int main()
{
	Matrix matrix;

	//Chapter2.1 估计5-20阶Hilbert矩阵的无穷范数；
	cout << endl;
	cout << "Chapter2.1 估计5-20阶Hilbert矩阵的无穷范数条件数" << endl;
	cout << endl;
	for (int n = 5; n <= 20; n++)
	{
		Hilbert_matrix(n, matrix); //生成n阶Hilbert矩阵；
		cout << n << "阶无穷范数条件数：" << Conditional_number_infinity(matrix) << endl;
		matrix.clear();
	}
	cout << endl;

	//Chapter2.2 随机选取x,计算b=Ax,然后用列主元Guass消去法求解方程组，估算5-30阶解的精度并与真实相对误差作比较；
	//内置函数选取的x伪随机数，其中选取x，默认为2^(1/2)*i，其次更换注释后，可以取2^(1/2)rand()；
	cout << endl;
	cout << "Chapter2.2 随机选取x,计算b=Ax,然后用列主元Guass消去法求解方程组，估算5-30阶解的精度并与真实相对误差作比较" << endl;
	cout << endl;
	for (int n = 5; n <= 30; n++)
	{
		vector<double> x(n), b(n), value(n),y(n);
		Vector_creat_random(n, x);                 //生成随机的x；
		creat_matrix(n, matrix);                   //生成系数矩阵A；
		cout << n << "阶所随机选取的x为：" << endl;
		Vector_cout(x);
		cout << endl;
		Matrix_times_vector(matrix, x, b);         //计算b；
		Guass_Line_solution(matrix, b, value);     //计算方程的解；
		Matrix_times_vector(matrix, value, y);     //计算A乘以计算解x；
		Vector_subtraction(y, b);                  //计算b-Ax并储存在y里面；
		double accuracy = Conditional_number_infinity(matrix) * vector_infinity_norm(y) / vector_infinity_norm(b);
		cout << n << "阶解的精度为：" << accuracy << endl;
		Vector_subtraction(value,x);               //计算精确解x与计算解的差值；
		double real_accuracy = vector_infinity_norm(value) / vector_infinity_norm(x);
		cout << n << "阶实际误差为：" << real_accuracy << endl << endl << endl;
		matrix.clear();
	}
}