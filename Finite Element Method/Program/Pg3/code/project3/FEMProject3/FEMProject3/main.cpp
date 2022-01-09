#include"FEMori.h"
#include"shinshkin.h"

using namespace std;

int main()
{
	int n = 7;
	vector<vector<double>> error(n, vector<double>(6, 0.0));
	error_out(error);
	/*
	cout << "基于一次有限元空间求解误差与收敛阶：" << endl;
	cout << "N   \tL无穷误差  \t误差收敛阶  \tL1误差  \t误差收敛阶" << endl;
	for (int i = 0; i < size(error); i++)
		cout << pow(2, i) * 10 << "\t" << error[i][0] << "\t" << error[i][1] << "\t" << error[i][2] << "\t" << error[i][3] << endl;
	cout << endl;
	*/
	
	cout << "基于shinshkin划分的一次有限元空间求解误差与收敛阶：" << endl;
	cout << "N\t time\t L无穷误差\t 误差收敛阶\t L1误差\t 误差收敛阶" << endl;
	shinshkin_error_out(n, error);
	matrix_output(error);

}