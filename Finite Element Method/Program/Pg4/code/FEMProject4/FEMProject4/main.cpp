#include"femori.h"
#include"ConjugateGradient.h"

using namespace std;

int main()
{
	for (int i = 0; i < 5; i++)
	{
		//导入数据；
		vector<vector<double>> node; vector<int> node_boundary_marker;
		vector<vector<int>> triangle; vector<vector<double>> tri_circumcenter;
		vector<vector<int>> edge;
		int file_i = i;
		data_input(file_i, node, node_boundary_marker, triangle, tri_circumcenter, edge);

		int n = node.size();;
		vector<vector<double>> cof_k(n, vector<double>(n, 0));
		vector<double> cof_f(n, 0), cof_u(n, 0);
		make_cof_k_f(cof_k, cof_f, node_boundary_marker, node, triangle, tri_circumcenter, edge);

		//共轭梯度法求解；
		int times = 0;             //共轭梯度迭代法迭代次数；
		double iteration_error = 0;//共轭梯度迭代法的输出参考误差；
		Conjugate_Gradient(n, cof_k, cof_f, 1e-20, cof_u, times, iteration_error);

		/*输出检查
		if (i == 0)
		{
			matrix_output(cof_k);
			vector_output(cof_f);
			vector_output(cof_u);
		}
		//*/

		error_max_2(file_i, cof_u, triangle, node, node_boundary_marker);
	}
	/*
	vector<vector<double>> matrix(15, vector<double>(15, 0));
	vector<double> x(15, 0);
	x[0] = 1;
	vector<double> f(15, 0);
	matrix[5][5] = 1; matrix[5][8] = 1; matrix[5][10] = 1;
	matrix[8][5] = 1; matrix[10][5] = 1;
	matrix[6][6] = 1; matrix[6][8] = 1; matrix[6][11] = 1;
	matrix[8][6] = 1; matrix[11][6] = 1;
	matrix[8][8] = 1; matrix[8][10] = 1; matrix[8][11] = 1;
	matrix[10][8] = 1; matrix[11][8] = 1;
	matrix[10][10] = 1; 
	matrix[11][11] = 1;
	f[5] = 3; f[6] = 3; f[8] = 5; f[10] = 3; f[11] = 3;
	for (int i = 0; i < 15; i++)
	{
		for (int j = 0; j < 15; j++)
			cout << matrix[i][j] << ",";
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < size(x); i++)
		cout << f[i] << ",";
	cout << endl << endl;
	int times = 0;             //共轭梯度迭代法迭代次数；
	double iteration_error = 0;//共轭梯度迭代法的输出参考误差；
	Conjugate_Gradient(15, matrix, f, 1e-20, x, times, iteration_error);

	cout << times << "    " << iteration_error << endl;
	for (int i = 0; i < size(x); i++)
		cout << x[i] << ",";
	*/
	/*
	cout << "基于一次有限元空间求解误差与收敛阶：" << endl;
	cout << "N   \tL无穷误差  \t误差收敛阶  \tL1误差  \t误差收敛阶" << endl;
	for (int i = 0; i < size(error); i++)
		cout << pow(2, i) * 10 << "\t" << error[i][0] << "\t" << error[i][1] << "\t" << error[i][2] << "\t" << error[i][3] << endl;
	cout << endl;
	/*
	cout << "基于二次有限元空间求解误差与收敛阶：" << endl;
	cout << "N\tL无穷误差\t误差收敛阶\tL1误差  \t误差收敛阶" << endl;
	for (int i = 0; i < size(error); i++)
		cout << pow(2, i) * 10 << "\t" << error[i][4] << "\t" << error[i][5] << endl;// << "\t" << error[i][6] << "\t" << error[i][7] << endl;

	cout << endl << exp(1) << " " << exp(2) << endl;
	cout << endl << log(exp(1)) << " " << log(exp(2)) << endl;
	*/

}