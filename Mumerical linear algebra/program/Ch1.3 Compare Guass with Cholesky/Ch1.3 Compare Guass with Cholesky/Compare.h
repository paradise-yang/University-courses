#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

int band_matrix(int n, vector<vector<double>>& matrix)      //生成带状系数方阵
{
	int i, j;
	vector<double> row(n);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i == j) row[j] = 10;
			else if (j == i + 1 || j == i - 1) row[j] = 1;
			else row[j] = 0;
		}
		matrix.push_back(row);
	}
	return (0);
}

int Hilbert_matrix(int n, vector<vector<double>>& matrix)  //生成Hilbert方阵；
{
	int i = 0, j = 0;
	vector<double> row(n);
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			row[j] = 1.0 / (i + j + 1);
		matrix.push_back(row);
	}
	return 0;
}

int lower_unit_triangular(int n, vector<vector<double>> matrix, vector<double>& value)  //单位下三角方阵前代法求解
{
	int i, j;
	for (i = 0; i < n - 1; i++)
		for (j = i + 1; j < n; j++)
			value[j] = value[j] - value[i] * matrix[j][i];
	return 0;
}

int lower_triangular(int n, vector<vector<double>> matrix, vector<double>& value)  //下三角方阵前代法求解
{
	int i, j;
	for (i = 0; i < n - 1; i++)
	{
		value[i] = value[i] / matrix[i][i];
		for (j = i + 1; j < n; j++)
			value[j] = value[j] - value[i] * matrix[j][i];
	}
	value[n - 1] = value[n - 1] / matrix[n - 1][n - 1];
	return 0;
}

int upper_unit_triangular(int n, vector<vector<double>> matrix, vector<double>& value)  //单位上三角方阵回代法求解
{
	int i, j;
	for (i = n - 1; i > 0; i--)
		for (j = 0; j < i; j++)
			value[j] = value[j] - value[i] * matrix[j][i];
	return 0;
}

int upper_triangular(int n, vector<vector<double>> matrix, vector<double>& value)  //上三角方阵回代法求解
{
	int i, j;
	for (i = n - 1; i > 0; i--)
	{
		value[i] = value[i] / matrix[i][i];
		for (j = 0; j < i; j++)
			value[j] = value[j] - value[i] * matrix[j][i];
	}
	value[0] = value[0] / matrix[0][0];
	return 0;
}

double Error_analysis(int n, vector<vector<double>>& matrix, vector<double>& value, vector<double>& exact)//误差分析；
{
	double error = 0, temp = 0;
	int i, j;
	for (i = 0; i < n; i++)
	{
		temp = 0;
		for (j = 0; j < n; j++)
			temp += matrix[i][j] * value[j];
		error += pow(exact[i] - temp, 2);
	}
	error = sqrt(error);
	return error;
}

int Guass_LU(int n, vector<vector<double>>& matrix)  //未选主元的Guass消去求解LU分解
{
	int k, i, j;
	for (k = 0; k < n - 1; k++)
	{
		for (i = k + 1; i < n; i++)
			matrix[i][k] = matrix[i][k] / matrix[k][k];
		for (i = k + 1; i < n; i++)
			for (j = k + 1; j < n; j++)
				matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
	}
	return 0;
}

int Guass_Line(int n, vector<vector<double>>& matrix, vector<double>& value)  //列主元Guass消去求解LU，并同时置换
{
	int k, i, j, p = 0;
	double max, temp;
	for (k = 0; k < n - 1; k++)
	{
		max = matrix[k][k];//寻找列主元
		p = k;
		for (i = k + 1; i < n; i++)
		{
			if (fabs(matrix[i][k]) > fabs(max))
			{
				max = matrix[i][k];
				p = i;
			}
		}
		if (p != k) //置换；
		{
			vector<double> pc(matrix[k]);
			matrix[k] = matrix[p];
			matrix[p] = pc;
			temp = value[k];
			value[k] = value[p];
			value[p] = temp;
		}

		if (abs(matrix[k][k]) !=0)//>= 1e-15)
		{
			for (i = k + 1; i < n; i++)
				matrix[i][k] = matrix[i][k] / matrix[k][k];
			for (j = k + 1; j < n; j++)
				for (i = k + 1; i < n; i++)
					matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
		}
	}
	return 0;
}