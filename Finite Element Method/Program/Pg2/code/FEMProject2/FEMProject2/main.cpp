#include"FEMori.h"

using namespace std;

int main()
{
	vector<vector<double>> error(5, vector<double>(8, 0.0));
	error_out(error);
	cout << "基于一次有限元空间求解误差与收敛阶：" << endl;
	cout << "N   \tL无穷误差  \t误差收敛阶  \tL1误差  \t误差收敛阶" << endl;
	for (int i = 0; i < size(error); i++)
		cout << pow(2, i) * 10 << "\t" << error[i][0] << "\t" << error[i][1] << "\t" << error[i][2] << "\t" << error[i][3] << endl;
	cout << endl;
	cout << "基于二次有限元空间求解误差与收敛阶：" << endl;
	cout << "N\tL无穷误差\t误差收敛阶\tL1误差  \t误差收敛阶" << endl;
	for (int i = 0; i < size(error); i++)
		cout << pow(2, i) * 10 << "\t" << error[i][4] << "\t" << error[i][5] << "\t" << error[i][6] << "\t" << error[i][7] << endl;

}