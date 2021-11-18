#include "Guass_GuassLine.h"
#include <ctime>

#define dim 84

using namespace std;

typedef double elem;

int main()
{
	int i, j, k;
	vector<vector<elem>> matrix,matrix2,matrix0;
	vector<elem> value(dim), exact(dim);
	double error = 0;
	clock_t start_time, end_time, start_time1, end_time1;

	band_matrix(dim, matrix0);

	band_matrix(dim, matrix);            //������ɴ�״ϵ������
	value[0] = 7;                        //��ʼ��ֵ��
	value[dim - 1] = 14; 
	for (i = 1; i < dim-1; ++i)
		value[i] = 15;
	exact = value;
	/*
	cout << endl << endl << endl;
	vector<vector<double>> a= { {0.05, 0.07, 0.06, 0.05} , {0.07, 0.10, 0.08, 0.07}, {0.06, 0.08, 0.10, 0.09}, {0.05, 0.07, 0.09, 0.10} };
	vector<double> b= { 0.23, 0.32, 0.33, 0.31 };
	vector<double> x = { 0,0,0,0 };
	Guass_LU(4, a);
	lower_triangular(4, a, b);
	upper_triangular(4, a, b);
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << a[i][j] << "\t";
		cout << endl;
	}
	for (j = 0; j < 4; j++)
		cout << b[j] << "  ";
	cout << endl << endl << endl;
	*/

	start_time = clock();                //��ȡ��һ��ʱ�䣻
	Guass_LU(dim, matrix);               //δѡ��Ԫ��Guass��ȥ��
	lower_triangular(dim, matrix, value);//�������ǵó�y;
	upper_triangular(dim, matrix, value);//�������ǵó�x;
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "δѡ��Ԫ�ĸ�˹��ȥ������Ϊ" << endl;
	for (j = 0; j < dim; j++)
		cout << value[j] << "  ";
	error = Error_analysis(dim, matrix0, value, exact);//��������
	cout << endl << endl << "������|Ax-b|=" << error << endl;
	cout << "The run time is " << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

	cout << endl << endl;
	
	band_matrix(dim, matrix2);            //�ٴ�������ɴ�״ϵ������
	value = exact;                        //�ٴγ�ʼ��ֵ��

	start_time1 = clock();//��ȡ��һ��ʱ�䣻
	Guass_Line(dim, matrix2, value);      //����ԪGuass��ȥ��
	lower_triangular(dim, matrix2, value);//�������ǵó�y;
	upper_triangular(dim, matrix2, value);//�������ǵó�x;
	end_time1 = clock();//��ȡ�ڶ���ʱ�䣻
	cout << "����Ԫ�ĸ�˹��ȥ������Ϊ" << endl;
	for (j = 0; j < dim; j++)
		cout << value[j] << " ";
	error = Error_analysis(dim, matrix0, value, exact);//��������
	cout << endl << endl << "������|Ax-b|=" << error << endl;
	cout << "The run time is " << (double)(end_time1 - start_time1) / CLOCKS_PER_SEC << endl;

	cout << endl << endl;

	cout << "test" << endl;
	vector<vector<double>> tempmatrix;
	Hilbert_matrix(7, tempmatrix);
	Guass_LU(7, tempmatrix);
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 7; j++)
			cout << tempmatrix[i][j] << "\t";
		cout << endl;
	}
	cout << endl << endl;

	vector<vector<double>> tempmatrix2;
	vector<double> tempvector(7, 1);
	Hilbert_matrix(7, tempmatrix2);
	Guass_Line(7, tempmatrix2,tempvector);
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 7; j++)
			cout << tempmatrix2[i][j] << "\t";
		cout << endl;
	}


	return 0;
}