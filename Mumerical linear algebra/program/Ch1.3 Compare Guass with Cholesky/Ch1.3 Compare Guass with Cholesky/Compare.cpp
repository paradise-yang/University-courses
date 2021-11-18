#include "Compare.h"
#include <ctime>

#define dim 100
#define dimh 40

using namespace std;

typedef double elem;

//对于括号一，选取的b与Cholesky分解一样，但需要取消注释！！！

int main()
{
	int i, j;
	vector<vector<elem>> matrix, matrix2, matrix0;
	vector<vector<elem>> Hilbert, Hilbert2, Hilbert0;
	vector<elem> value(dim), exact(dim);
	double error = 0;
	clock_t start_time, end_time;

	cout << endl << "First part:随机选取b，用未选主元Guass消去法与列主元Guass消去法求解方程组。" << endl << endl;

	band_matrix(dim, matrix0);           //记录系数矩阵，便于误差分析；

	band_matrix(dim, matrix);            //生成题干带状系数方阵；
	/*for (i = 0; i < dim; i++)          //随机初始化值；
		value[i] = rand();*/
	value[0] = 11;                       //特殊“1”初始化值；
	value[dim - 1] = 11;
	for (i = 1; i < dim - 1; ++i)
		value[i] = 12;
	exact = value;                       //记录初始化，便于误差分析；
	cout << "随机选取的b为：" << endl;
	for (i = 0; i < dim; i++)
		cout << exact[i] << "  ";
	cout << endl << endl;

	start_time = clock();                //获取第一个时间；
	Guass_LU(dim, matrix);                    //未选主元Guass消去计算LU分解；
	lower_unit_triangular(dim, matrix, value);//解单位下三角得出y;
	upper_triangular(dim, matrix, value);     //解上三角得出x;
	end_time = clock();                  //获取第二个时间；
	cout << "未选主元Guass消去法，解为" << endl;
	for (j = 0; j < dim; j++)
		cout << value[j] << " ";
	error = Error_analysis(dim, matrix0, value, exact);//误差分析；
	cout << endl << endl << "误差分析|Ax-b|=" << error << endl;
	cout << "The run time is " << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

	cout << endl << endl;

	band_matrix(dim, matrix2);                 //再次生成题干带状系数方阵；
	value = exact;                             //再次初始化值；

	start_time = clock();                //获取第一个时间；
	Guass_Line(dim, matrix2,value);            //列主元Guass消去计算LU分解；
	lower_unit_triangular(dim, matrix2, value);//解单位下三角得出y;
	upper_triangular(dim, matrix2, value);     //解上三角得出x;
	end_time = clock();                  //获取第二个时间；
	cout << "列主元Guass消去法，解为" << endl;
	for (j = 0; j < dim; j++)
		cout << value[j] << " ";
	error = Error_analysis(dim, matrix0, value, exact);//误差分析；
	cout << endl << endl << "误差分析|Ax-b|=" << error << endl;
	cout << "The run time is " << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

	cout << endl << endl;

	cout << endl << "Second part:构造Hilbert方阵，用未选主元Guass消去法与列主元Guass消去法求解Hilbert方程组。" << endl << endl;

	Hilbert_matrix(dimh, Hilbert0);           //记录系数矩阵，便于误差分析；

	Hilbert_matrix(dimh, Hilbert);            //生成题干Hilbert系数方阵；
	for (i = 0; i < dimh; i++)                //初始化值；
	{
		value[i] = 0;
		for (j = 0; j < dimh; j++)
			value[i] += Hilbert[i][j];
	}
	exact = value;                         //记录初始化，便于误差分析；

	start_time = clock();                //获取第一个时间；
	Guass_LU(dimh, Hilbert);                    //未选主元Guass消去法计算LU分解；
	lower_unit_triangular(dimh, Hilbert, value);//解单位下三角得出y;
	upper_triangular(dimh, Hilbert, value);     //解上三角得出x;
	end_time = clock();                  //获取第二个时间；
	cout << "未选主元Guass消去法，解为" << endl;
	for (j = 0; j < dimh; j++)
		cout << value[j] << " ";
	error = Error_analysis(dimh, Hilbert0, value, exact);//误差分析；
	cout << endl << endl << "误差分析|Ax-b|=" << error << endl;
	cout << "The run time is " << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

	cout << endl << endl;

	Hilbert_matrix(dimh, Hilbert2);              //再次生成题干Hilbert系数方阵；
	value = exact;                               //再次初始化值；

	start_time = clock();                //获取第一个时间；
	Guass_Line(dimh, Hilbert2,value);            //列主元Guass消去法计算LU分解；
	lower_unit_triangular(dimh, Hilbert2, value);//解单位下三角得出y;
	upper_triangular(dimh, Hilbert2, value);     //解上三角得出x;
	end_time = clock();                  //获取第二个时间；
	cout << "列主元Guass消去法，解为" << endl;
	for (j = 0; j < dimh; j++)
		cout << value[j] << " ";
	error = Error_analysis(dimh, Hilbert0, value, exact);//误差分析；
	cout << endl << endl << "误差分析|Ax-b|=" << error << endl;
	cout << "The run time is " << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

	cout << endl << endl;

	return 0;
}

/*for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++)
			cout << matrix[i][j]<<" ";
		cout << endl;
	}输出矩阵函数*/