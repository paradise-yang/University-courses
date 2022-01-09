#include"uniform.h"
#include"shinshkin.h"

using namespace std;

int main()
{
	int n = 8;
	vector<vector<double>> error(n, vector<double>(8, 0.0));
	/* 
	cout << "\epsilon=0.1，均匀网格。" << endl << endl;

	cout << "基于一次有限元空间求解误差与收敛阶：" << endl;
	cout << "N\t time\t L无穷误差\t 误差收敛阶\t L1误差\t 误差收敛阶\t L2误差\t误差收敛阶" << endl;
	CG_error_out(error, n, epsilon1, 1);
	matrix_output(error);
	cout << "基于二次有限元空间求解误差与收敛阶：" << endl;
	cout << "N\t time\t L无穷误差\t 误差收敛阶\t L1误差\t 误差收敛阶\t L2误差\t误差收敛阶" << endl;
	CG_error_out_2(error, n, epsilon1, 1);
	matrix_output(error);
	cout << endl;
	cout << "*******************************************************************************************************************" << endl;
	cout << endl << endl;
	//*/
	/* 
	cout << "\epsilon=1e-7，均匀网格。" << endl << endl;

	cout << "基于一次有限元空间求解误差与收敛阶：" << endl;
	cout << "N\t time\t L无穷误差\t 误差收敛阶\t L1误差\t 误差收敛阶\t L2误差\t误差收敛阶" << endl;
	CG_error_out(error, n, epsilon7, 7);
	matrix_output(error);
	cout << "基于二次有限元空间求解误差与收敛阶：" << endl;
	cout << "N\t time\t L无穷误差\t 误差收敛阶\t L1误差\t 误差收敛阶\t L2误差\t误差收敛阶" << endl;
	CG_error_out_2(error, n, epsilon7, 7);
	matrix_output(error);
	cout << endl;
	cout << "*******************************************************************************************************************" << endl;
	cout << endl << endl;
	//*/
	///* 
	cout << "\epsilon=1e-7，shinshkin网格。" << endl << endl;

	cout << "基于shinshkin划分的一次有限元空间求解误差与收敛阶：" << endl;
	cout << "N\t time\t L无穷误差\t 误差收敛阶\t L1误差\t 误差收敛阶" << endl;
	shinshkin_error_out(n, error, 7);
	matrix_output(error);
	cout << "基于shinshkin划分的二次有限元空间求解误差与收敛阶：" << endl;
	cout << "N\t time\t L无穷误差\t 误差收敛阶\t L1误差\t 误差收敛阶" << endl;
	shinshkin_error_out_2(n, error, 7);
	matrix_output(error);
	cout << endl;
	//*/
}