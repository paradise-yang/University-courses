#include"Cholesky_CholeskyImprove.h"
#include <ctime>

#define dim 100
#define dimh 10

using namespace std;

typedef double elem;

//对于第二题第一问，给出两种b的选取，一种是使得解均为1的b，一种是用系统内置函数生成随机数，但需要取消注释运行额外一种。

int main()
{
	int i, j;
	vector<vector<elem>> matrix, matrix2, matrix0;
	vector<vector<elem>> Hilbert, Hilbert2, Hilbert0;
	vector<elem> value(dim), exact(dim);
	double error = 0;
	clock_t start_time, end_time, start_time1, end_time1;

	cout << endl << "First part:随机选取b，用平方根法与改进平方根法求解方程组。" << endl << endl;
	
	band_matrix(dim, matrix0);           //记录系数矩阵，便于误差分析；

	band_matrix(dim, matrix);            //生成题干带状系数方阵；
	/*for (i = 0; i < dim; i++)            //随机初始化值；
		value[i] = rand();*/
	value[0] = 11;                         //特殊“1”初始化值；
	value[dim - 1] = 11;
	for (i = 1; i < dim-1; ++i)
		value[i] = 12;
	exact = value;                       //记录初始化，便于误差分析；
	cout << "随机选取的b为：" << endl;
	for (i = 0; i < dim; i++)
		cout << exact[i] << "  ";
	cout << endl << endl;
	
	start_time = clock();                //获取第一个时间；
	Cholesky(dim, matrix);               //平方根法计算Cholesky分解；
	lower_triangular(dim, matrix, value);//解下三角得出y;
	upper_triangular(dim, matrix, value);//解上三角得出x;
	end_time = clock();                  //获取第二个时间；
	cout << "平方根法Cholesky分解，解为" << endl;
	for (j = 0; j < dim; j++)
		cout << value[j] << " ";
	error = Error_analysis(dim, matrix0, value, exact);//误差分析；
	cout << endl << endl << "误差分析|Ax-b|=" << error << endl;
	cout << "The run time is " << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

	cout << endl << endl;

	band_matrix(dim, matrix2);                 //再次生成题干带状系数方阵；
	value = exact;                             //再次初始化值；

	start_time = clock();                //获取第一个时间；
	Cholesky_improve(dim, matrix2);            //改进平方根法计算Cholesky分解；
	lower_unit_triangular(dim, matrix2, value);//解单位下三角得出y;
	D(dim, matrix2, value);                    //计算对角阵得出z;
	upper_unit_triangular(dim, matrix2, value);//解单位上三角得出x;
	end_time = clock();                  //获取第二个时间；
	cout << "改进平方根法Cholesky分解，解为" << endl;
	for (j = 0; j < dim; j++)
		cout << value[j] << " ";
	error = Error_analysis(dim, matrix0, value, exact);//误差分析；
	cout << endl << endl << "误差分析|Ax-b|=" << error << endl;
	cout << "The run time is " << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

	cout << endl << endl;

	cout << endl << "Second part:构造Hilbert方阵，用平方根法与改进平方根法求解方程组。" << endl << endl;

	Hilbert_matrix(dimh, Hilbert0);           //记录系数矩阵，便于误差分析；

	Hilbert_matrix(dimh, Hilbert);            //生成题干Hilbert系数方阵；
	for (i = 0; i < dimh; i++)                //初始化值；
	{
		value[i] = 0;
		for (j = 0; j < dimh; j++)
			value[i] += Hilbert[i][j];
	}
	exact.clear();
	exact = value;                         //记录初始化，便于误差分析；

	start_time1 = clock();//获取第一个时间；
	Cholesky(dimh, Hilbert);               //平方根法计算Cholesky分解；
	lower_triangular(dimh, Hilbert, value);//解下三角得出y;
	upper_triangular(dimh, Hilbert, value);//解上三角得出x;
	end_time1 = clock();//获取第二个时间；
	cout << "平方根法Cholesky分解，解为" << endl;
	for (j = 0; j < dimh; j++)
		cout << value[j] << " ";
	error = Error_analysis(dimh, Hilbert0, value, exact);//误差分析；
	cout << endl << endl << "误差分析|Ax-b|=" << error << endl;
	cout << "The run time is " << (double)(end_time1 - start_time1) / CLOCKS_PER_SEC << endl;

	cout << endl << endl;

	Hilbert_matrix(dimh, Hilbert2);                 //再次生成题干Hilbert系数方阵；
	value.clear();
	value = exact;                               //再次初始化值；
	//for (j = 0; j < dimh; j++)
	//	cout << value[j] << " ";
	//cout << endl << endl;

	start_time1 = clock();//获取第一个时间；
	Cholesky_improve(dimh, Hilbert2);            //改进平方根法计算Cholesky分解；
	lower_unit_triangular(dimh, Hilbert2, value); //解单位下三角得出y;
	D(dimh, Hilbert2, value);                    //计算对角阵得出z;
	upper_unit_triangular(dimh, Hilbert2, value);//解单位上三角得出x;
	end_time1 = clock();//获取第二个时间；
	cout << "改进平方根法Cholesky分解，解为" << endl;
	for (j = 0; j < dimh; j++)
		cout << value[j] << " ";
	error = Error_analysis(dimh, Hilbert0, value, exact);//误差分析；
	cout << endl << endl << "误差分析|Ax-b|=" << error << endl;
	cout << "The run time is " << (double)(end_time1 - start_time1) / CLOCKS_PER_SEC << endl;
	
	cout << endl << endl;

	return 0;
}

/*for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++)
			cout << matrix2[i][j]<<" ";
		cout << endl;
	}
	for (j = 0; j < dim; j++)
		cout << value[j] << " ";
	cout << endl << endl;*/