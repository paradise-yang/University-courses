#include "Compare.h"
#include <ctime>

#define dim 100
#define dimh 40

using namespace std;

typedef double elem;

//��������һ��ѡȡ��b��Cholesky�ֽ�һ��������Ҫȡ��ע�ͣ�����

int main()
{
	int i, j;
	vector<vector<elem>> matrix, matrix2, matrix0;
	vector<vector<elem>> Hilbert, Hilbert2, Hilbert0;
	vector<elem> value(dim), exact(dim);
	double error = 0;
	clock_t start_time, end_time;

	cout << endl << "First part:���ѡȡb����δѡ��ԪGuass��ȥ��������ԪGuass��ȥ����ⷽ���顣" << endl << endl;

	band_matrix(dim, matrix0);           //��¼ϵ�����󣬱�����������

	band_matrix(dim, matrix);            //������ɴ�״ϵ������
	/*for (i = 0; i < dim; i++)          //�����ʼ��ֵ��
		value[i] = rand();*/
	value[0] = 11;                       //���⡰1����ʼ��ֵ��
	value[dim - 1] = 11;
	for (i = 1; i < dim - 1; ++i)
		value[i] = 12;
	exact = value;                       //��¼��ʼ����������������
	cout << "���ѡȡ��bΪ��" << endl;
	for (i = 0; i < dim; i++)
		cout << exact[i] << "  ";
	cout << endl << endl;

	start_time = clock();                //��ȡ��һ��ʱ�䣻
	Guass_LU(dim, matrix);                    //δѡ��ԪGuass��ȥ����LU�ֽ⣻
	lower_unit_triangular(dim, matrix, value);//�ⵥλ�����ǵó�y;
	upper_triangular(dim, matrix, value);     //�������ǵó�x;
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "δѡ��ԪGuass��ȥ������Ϊ" << endl;
	for (j = 0; j < dim; j++)
		cout << value[j] << " ";
	error = Error_analysis(dim, matrix0, value, exact);//��������
	cout << endl << endl << "������|Ax-b|=" << error << endl;
	cout << "The run time is " << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

	cout << endl << endl;

	band_matrix(dim, matrix2);                 //�ٴ�������ɴ�״ϵ������
	value = exact;                             //�ٴγ�ʼ��ֵ��

	start_time = clock();                //��ȡ��һ��ʱ�䣻
	Guass_Line(dim, matrix2,value);            //����ԪGuass��ȥ����LU�ֽ⣻
	lower_unit_triangular(dim, matrix2, value);//�ⵥλ�����ǵó�y;
	upper_triangular(dim, matrix2, value);     //�������ǵó�x;
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "����ԪGuass��ȥ������Ϊ" << endl;
	for (j = 0; j < dim; j++)
		cout << value[j] << " ";
	error = Error_analysis(dim, matrix0, value, exact);//��������
	cout << endl << endl << "������|Ax-b|=" << error << endl;
	cout << "The run time is " << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

	cout << endl << endl;

	cout << endl << "Second part:����Hilbert������δѡ��ԪGuass��ȥ��������ԪGuass��ȥ�����Hilbert�����顣" << endl << endl;

	Hilbert_matrix(dimh, Hilbert0);           //��¼ϵ�����󣬱�����������

	Hilbert_matrix(dimh, Hilbert);            //�������Hilbertϵ������
	for (i = 0; i < dimh; i++)                //��ʼ��ֵ��
	{
		value[i] = 0;
		for (j = 0; j < dimh; j++)
			value[i] += Hilbert[i][j];
	}
	exact = value;                         //��¼��ʼ����������������

	start_time = clock();                //��ȡ��һ��ʱ�䣻
	Guass_LU(dimh, Hilbert);                    //δѡ��ԪGuass��ȥ������LU�ֽ⣻
	lower_unit_triangular(dimh, Hilbert, value);//�ⵥλ�����ǵó�y;
	upper_triangular(dimh, Hilbert, value);     //�������ǵó�x;
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "δѡ��ԪGuass��ȥ������Ϊ" << endl;
	for (j = 0; j < dimh; j++)
		cout << value[j] << " ";
	error = Error_analysis(dimh, Hilbert0, value, exact);//��������
	cout << endl << endl << "������|Ax-b|=" << error << endl;
	cout << "The run time is " << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

	cout << endl << endl;

	Hilbert_matrix(dimh, Hilbert2);              //�ٴ��������Hilbertϵ������
	value = exact;                               //�ٴγ�ʼ��ֵ��

	start_time = clock();                //��ȡ��һ��ʱ�䣻
	Guass_Line(dimh, Hilbert2,value);            //����ԪGuass��ȥ������LU�ֽ⣻
	lower_unit_triangular(dimh, Hilbert2, value);//�ⵥλ�����ǵó�y;
	upper_triangular(dimh, Hilbert2, value);     //�������ǵó�x;
	end_time = clock();                  //��ȡ�ڶ���ʱ�䣻
	cout << "����ԪGuass��ȥ������Ϊ" << endl;
	for (j = 0; j < dimh; j++)
		cout << value[j] << " ";
	error = Error_analysis(dimh, Hilbert0, value, exact);//��������
	cout << endl << endl << "������|Ax-b|=" << error << endl;
	cout << "The run time is " << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

	cout << endl << endl;

	return 0;
}

/*for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++)
			cout << matrix[i][j]<<" ";
		cout << endl;
	}���������*/