#include "QR_LS.h"

#define dim1 30   //Ch1.1的维数，验证发现30阶之后开始有一定偏差，50阶起误差偏大，55是临界值，当维数 >55 时只输出 -nan(ind)；
#define dim2 100  //Ch1.2.(1)的维数；
#define dim3 10   //Ch1.2.(2)的维数，验证发现11是临界值，当维数 >11 时输出结果已经开始偏差较大；

using namespace std;

typedef  vector<vector<double>> Matrix;

//特注：给出了两种算法，
//一种是老师所说的“偷懒”生成实际Householder阵运用矩阵乘法计算QR分解与LS问题，
//另一种是不直接生成Householder阵通过其生成向量与系数计算QR分解，
//这里由于注释原因，只采用了第一种解法，如果想看第二种运算结果，需在头文件中更改注释！！！

int main()
{
	Matrix matrix, matrix2;

	//Ch3.1 QR方法求解Ch1的三个方程组，并将所得结果进行比较（比较在实验报告里面给出，这里不再重新输出Ch1的计算结果）；
	cout << endl << "Ch3.1 QR方法求解Ch1的三个方程组，并将所得结果进行比较." << endl;

	band_matrix1(dim1, matrix);//初始化矩阵；
	vector<double> b1(dim1);
	b1[0] = 7;                 //初始化向量值；
	b1[dim1-1] = 14;
	for (int i = 1; i < dim1-1; ++i)
		b1[i] = 15;
	LS(matrix, dim1, dim1, b1);
	cout << endl << "Ch1.1的线性方程组QR解得结果为：" << endl;
	Vector_cout(dim1,b1);
	matrix.clear();

	band_matrix2(dim2, matrix);//初始化矩阵；
	vector<double> b2(dim2);
	b2[0] = 11;                //初始化向量值；
	b2[dim2-1] = 11;
	for (int i = 1; i < dim2-1; ++i)
		b2[i] = 12;
	LS(matrix, dim2, dim2, b2);
	cout << endl << "Ch1.2.(1)的线性方程组QR解得结果为：" << endl;
	Vector_cout(dim2,b2);
	matrix.clear();

	Hilbert_matrix(dim3, matrix);
	vector<double> b3(dim3);
	for (int i = 0; i < dim3; i++)  //初始化向量值；
	{
		b3[i] = 0;
		for (int j = 0; j < dim3; j++)
			b3[i] += matrix[i][j];
	}
	LS(matrix, dim3, dim3, b3);
	cout << endl << "Ch1.2.(2)的线性方程组QR解得结果为：" << endl;
	Vector_cout(dim3,b3);
	cout << endl;
	matrix.clear();


	//Ch3.2 求二次多项式，使得残向量在的2-范数最小意义下拟合相应数据；
	cout << endl << endl << "Ch3.2 求二次多项式，使得残量在2-范数最小意义下拟合相应数据." << endl;
	vector<double> b4(7), y4(7);
	Ch3_2_make(matrix, b4);
	matrix2 = matrix; y4 = b4;
	LS(matrix, 7, 3, b4);
	cout << endl << "所求多项式为：y=" << b4[0] << " t^2 + " << b4[1] << " t + " << b4[2];
	cout << endl << "此时残向量2-范数为：||b - Ax|| = " << sqrt(norm_2(matrix2, 7, 3, y4, b4)) << endl;
	cout << endl;
	matrix.clear(); 
	matrix2.clear();


	//Ch3.3 给出房产估计的线性模型最小二乘解；
	cout << endl << endl << "Ch3.3 求房产估计的线性模型最小二乘解." << endl;
	vector<double> b5(28), y5(28);
	Ch3_3_make(matrix, b5);
	matrix2 = matrix; y5 = b5;
	LS(matrix, 28, 12, b5);
	cout << endl << "该线性模型的最小二乘解为：" << endl;
	Vector_cout(12, b5);
	cout << "此时残向量2-范数为：||b - Ax|| = " << sqrt(norm_2(matrix2, 28, 12, y5, b5)) << endl;
	cout << endl;
	matrix.clear();
	matrix2.clear();
}