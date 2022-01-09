#include"femori.h"

using namespace std;

int main()
{
	int n = 9;
	vector<vector<double>> error(n, vector<double>(6, 0));
	/*
	cout << "共轭梯度法求解：" << endl;
	cout << "N\t运行时间\t L无穷误差\t误差收敛阶\t L1误差\t 误差收敛阶" << endl;
	CG_error_out(error, n);
	matrix_output(error);
	*/

	cout << "多重网格法求解：" << endl;
	cout << "N\t运行时间\t L无穷误差\t误差收敛阶\t L1误差\t误差收敛阶" << endl;
	vector<vector<double>> errors(n, vector<double>(6, 0));
	MG_error_out_two(errors, n);
	matrix_output(errors);

	/*
	cout << "多重网格two法求解：" << endl;
	cout << "N\t运行时间\t L无穷误差\t误差收敛阶\t L1误差\t误差收敛阶" << endl;
	vector<vector<double>> errorss(n, vector<double>(6, 0));
	MG_error_out_two(errorss, n);
	matrix_output(errorss);
	//*/
}