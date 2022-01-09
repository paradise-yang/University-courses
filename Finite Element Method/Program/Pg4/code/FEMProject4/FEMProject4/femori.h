#pragma once
 
#include<iostream>
#include <iomanip>
#include <vector>
#include<cmath>
#include <fstream>  //�ļ����⺯��
#include <string>

#define pi 3.141592653589793238462643383279502884197169399375105820974944592307816

using namespace std;

//�����i���ļ�����ض��㡢�����Ρ��߼����ݣ�
void data_input(int i, vector<vector<double>> &node, vector<int> &node_boundary_marker, vector<vector<int>> &triangle, vector<vector<double>> &tri_circumcenter, vector<vector<int>> &edge)
{
	vector<double> temp(1, 0);  //����
	ifstream infile;        //������
	ofstream outfile;       //�����
	string path = "E:\\study_materials\\FEM\\Pg4\\MESH\\gd";

	//����ڵ㣻
	infile.open(path + to_string(i) + ".n", ios::in);
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	int node_n = 0;//�ڵ����   
	infile >> node_n;
	vector<double> temp_node(2, 0);//�ڵ��������¼�ڵ����ꣻ
	for (int i = 0; i < node_n; i++)
		infile >> temp[0], infile >> temp_node[0], infile >> temp_node[1], infile >> temp[0],
		node_boundary_marker.push_back(temp[0]),
		node.push_back(temp_node);
	infile.close();   //�ر�data�ļ�
	/*������
	for (int i = 0; i < node.size(); i++)
		cout <<"{"<< node[i][0] << "," << node[i][1]<<"},";
	cout << endl << endl;
	//*/
	//���������Σ�
	infile.open(path + to_string(i) + ".e", ios::in);
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	int tri_n = 0;//�����θ�����
	infile >> tri_n;
	vector<int> temp_triangle(9, 0);//��������صĵ㡢�������Ρ��ߵ���Ϣ�����������Σ�λ��4��5��6�У���Ϊ-1��ò�Ϊ�߽��������Σ�
	vector<double> temp_tri_circumcenter(2, 0);//�������������ꣻ
	for (int i = 0; i < tri_n; i++)
	{
		infile >> temp[0];
		for (int j = 0; j < 9; j++)
			infile >> temp_triangle[j];
		triangle.push_back(temp_triangle);
		infile >> temp_tri_circumcenter[0], infile >> temp_tri_circumcenter[1];
		tri_circumcenter.push_back(temp_tri_circumcenter);
		infile >> temp[0];
	}
	infile.close();   //�ر�data�ļ�
	/*������
	for (int i = 0; i < triangle.size(); i++)
	{
		for (int j = 0; j < triangle[0].size(); j++)
			cout << triangle[i][j] << "\t";
		cout << tri_circumcenter[i][0] << "\t" << tri_circumcenter[i][1] << endl;
	}
	//*/

	//����ߣ�
	infile.open(path + to_string(i) + ".s", ios::in);
	if (!infile.is_open())
		cout << "Open file failure" << endl;
	int edge_n = 0;//�߸�����
	infile >> edge_n;
	vector<int> temp_edge(4, 0);//�ߵ���ʼ���սᶥ�㣬������������ţ�����4��ֵΪ-1���Ҳ��������Σ�
	for (int i = 0; i < edge_n; i++)
	{
		infile >> temp[0];
		for (int j = 0; j < 4; j++)
			infile >> temp_edge[j];
		edge.push_back(temp_edge);
		infile >> temp[0];
	}
	infile.close();   //�ر�data�ļ�
	/*������
	for (int i = 0; i < edge.size(); i++)
	{
		for (int j = 0; j < edge[0].size(); j++)
			cout << edge[i][j] << "\t";
		cout << endl;
	}
	//*/
}

//u(x,y);
double u(double x, double y)
{
	return sin(pi*x)*sin(pi*y);
}

//f(x,y);
double f(double x, double y)
{
	return 2 * pi*pi*sin(pi*x)*sin(pi*y);
}

//װ����������Ӧ���㣬�����Ե�һ����Ϊ���ģ����ɶ�Ӧa(v,v)ֵ��
double diagnoal_integrate(int node_i, int node_j, int node_k, vector<vector<double>> node)
{
	vector<vector<double>> temp(2, vector<double>(2, 0));
	temp[0][0] = node[node_j][0] - node[node_i][0], temp[0][1] = node[node_k][0] - node[node_i][0];
	temp[1][0] = node[node_j][1] - node[node_i][1], temp[1][1] = node[node_k][1] - node[node_i][1];
	double fabs_det = fabs(temp[0][0] * temp[1][1] - temp[0][1] * temp[1][0]);

	return ((temp[0][0] - temp[0][1])*(temp[0][0] - temp[0][1]) + (temp[1][0] - temp[1][1])*(temp[1][0] - temp[1][1])) / 2 / fabs_det;
}

//װ���������ָ�����㣬��ǰ������Ϊָ�������ɵĻ������������Ӧ��a(u,v)ֵ��
double internal_integrate(int node_i, int node_j, int node_k, vector<vector<double>> node)
{
	vector<vector<double>> temp(2, vector<double>(2, 0));
	temp[0][0] = node[node_j][0] - node[node_i][0], temp[0][1] = node[node_k][0] - node[node_i][0];
	temp[1][0] = node[node_j][1] - node[node_i][1], temp[1][1] = node[node_k][1] - node[node_i][1];
	double fabs_det = fabs(temp[0][0] * temp[1][1] - temp[0][1] * temp[1][0]);

	return ((temp[0][0] - temp[0][1])*temp[0][1] + (temp[1][0] - temp[1][1])*temp[1][1]) / 2 / fabs_det;
}

//ת�����������淶��Ԫ�µĵ�ת����һ�㵥Ԫ�µ����ꣻ
void transform_normal_to_common(int node_i, int node_j, int node_k, vector<vector<double>> node, double& x, double& y)
{
	vector<vector<double>> temp(2, vector<double>(2, 0));
	temp[0][0] = node[node_j][0] - node[node_i][0], temp[0][1] = node[node_k][0] - node[node_i][0];
	temp[1][0] = node[node_j][1] - node[node_i][1], temp[1][1] = node[node_k][1] - node[node_i][1];
	double temp_x = x, temp_y = y;
	x = node[node_i][0] + temp[0][0] * temp_x + temp[0][1] * temp_y;
	y = node[node_i][1] + temp[1][0] * temp_x + temp[1][1] * temp_y;
}

//װ������������ָ��ɢ���㣬�����Ե�һ���ڵ�Ϊ���ģ����ɶ�ӦF(v)ֵ�����λ���
double make_f_integrate(int node_i, int node_j, int node_k, vector<vector<double>> node)
{
	vector<vector<double>> temp(2, vector<double>(2, 0));
	temp[0][0] = node[node_j][0] - node[node_i][0], temp[0][1] = node[node_k][0] - node[node_i][0];
	temp[1][0] = node[node_j][1] - node[node_i][1], temp[1][1] = node[node_k][1] - node[node_i][1];
	double fabs_det = fabs(temp[0][0] * temp[1][1] - temp[0][1] * temp[1][0]);
	//double temp = 3 * f(node[node_i][0], node[node_i][1]);
	return (3 * f(node[node_i][0], node[node_i][1])
		+ 4 * f((node[node_i][0] + node[node_j][0]) / 2, (node[node_i][1] + node[node_j][1]) / 2)
		+ 4 * f((node[node_i][0] + node[node_k][0]) / 2, (node[node_i][1] + node[node_k][1]) / 2)
		+ 9 * f((node[node_i][0] + node[node_j][0] + node[node_k][0]) / 3, (node[node_i][1] + node[node_j][1] + node[node_k][1]) / 3))*fabs_det / 120;
	//return f(node[node_i][0], node[node_i][1])*fabs_det / 4;
}

//�������ε�Ԫɨ�裬װ��նȾ��󣬲�װ���Ҳ�ϵ��������
void make_cof_k_f(vector<vector<double>>& cof_k, vector<double>& cof_f, vector<int> node_boundary_marker, vector<vector<double>> node, vector<vector<int>> triangle, vector<vector<double>> tri_circumcenter, vector<vector<int>> edge)
{
	int n = node.size();//�նȾ����ά��Ϊ���������
	for (int e = 0; e < triangle.size(); e++)
	{
		if (!node_boundary_marker[triangle[e][0]])//�Ǳ߽�������
		{
			cof_k[triangle[e][0]][triangle[e][0]] += diagnoal_integrate(triangle[e][0], triangle[e][1], triangle[e][2], node);
			cof_f[triangle[e][0]] += make_f_integrate(triangle[e][0], triangle[e][1], triangle[e][2], node);
			if (!node_boundary_marker[triangle[e][1]])
				cof_k[triangle[e][0]][triangle[e][1]] += internal_integrate(triangle[e][0], triangle[e][1], triangle[e][2], node),
				cof_k[triangle[e][1]][triangle[e][0]] += internal_integrate(triangle[e][0], triangle[e][1], triangle[e][2], node);
			if (!node_boundary_marker[triangle[e][2]])
				cof_k[triangle[e][0]][triangle[e][2]] += internal_integrate(triangle[e][2], triangle[e][0], triangle[e][1], node),
				cof_k[triangle[e][2]][triangle[e][0]] += internal_integrate(triangle[e][2], triangle[e][0], triangle[e][1], node);
		}
		if (!node_boundary_marker[triangle[e][1]])
		{
			cof_k[triangle[e][1]][triangle[e][1]] += diagnoal_integrate(triangle[e][1], triangle[e][2], triangle[e][0], node);
			cof_f[triangle[e][1]] += make_f_integrate(triangle[e][1], triangle[e][2], triangle[e][0], node);
			if (!node_boundary_marker[triangle[e][2]])
				cof_k[triangle[e][1]][triangle[e][2]] += internal_integrate(triangle[e][1], triangle[e][2], triangle[e][0], node),
				cof_k[triangle[e][2]][triangle[e][1]] += internal_integrate(triangle[e][1], triangle[e][2], triangle[e][0], node);
		}
		if(!node_boundary_marker[triangle[e][2]])
			cof_k[triangle[e][2]][triangle[e][2]] += diagnoal_integrate(triangle[e][2], triangle[e][0], triangle[e][1], node),
			cof_f[triangle[e][2]] += make_f_integrate(triangle[e][2], triangle[e][0], triangle[e][1], node);

	}
}

//�������������������ֻ�����������true �� ��false��
bool out_product(int node_i, int node_j, vector<vector<double>> node, double x, double y)
{
	if ((x - node[node_i][0])*(node[node_j][1] - node[node_i][1]) - (node[node_j][0] - node[node_i][0])*(y - node[node_i][1]) > 0)
		return true;
	else
		return false;
}

//���������������ڵ�λ�ã��жϵ�(x,y)�Ƿ����������ڣ���true������false;
bool inside_triangle(int i, int j, int k, vector<vector<double>> node, double x, double y)
{
	bool first = out_product(i, j, node, x, y);
	bool second = out_product(j, k, node, x, y);
	bool third = out_product(k, i, node, x, y);
	if (first == second && second == third)
		return true;
	else 
		return false;
}

//��������������i�ڵ�Ϊ���ĺ����ĺ���ֵbasis(x,y);
double basis_i(int node_i, int node_j, int node_k, vector<vector<double>> node, double x, double y)
{
	vector<vector<double>> temp(2, vector<double>(2, 0));
	temp[0][0] = node[node_j][0] - node[node_i][0], temp[0][1] = node[node_k][0] - node[node_i][0];
	temp[1][0] = node[node_j][1] - node[node_i][1], temp[1][1] = node[node_k][1] - node[node_i][1];
	double det = temp[0][0] * temp[1][1] - temp[0][1] * temp[1][0];
	
	return ((x - node[node_i][0])*(node[node_j][1] - node[node_k][1]) + (y - node[node_i][1])*(node[node_k][0] - node[node_j][0])) / det + 1;
}


//ͨ��ϵ��������ֵ��⺯��uh(x,y);
double uh(double x, double y, vector<double> cof_u, vector<vector<int>> triangle, vector<vector<double>> node, vector<int> node_boundary_marker)
{
	double temp = 0;
	for (int e = 0; e < triangle.size(); e++)
	{
		if (inside_triangle(triangle[e][0], triangle[e][1], triangle[e][2], node, x, y))
		{//��(x,y)���������ڲ��в�ѯ�������Ҫ��
			if (!node_boundary_marker[triangle[e][0]])
				temp += cof_u[triangle[e][0]] * basis_i(triangle[e][0], triangle[e][1], triangle[e][2], node, x, y);
			if (!node_boundary_marker[triangle[e][1]])
				temp += cof_u[triangle[e][1]] * basis_i(triangle[e][1], triangle[e][2], triangle[e][0], node, x, y);
			if (!node_boundary_marker[triangle[e][2]])
				temp += cof_u[triangle[e][2]] * basis_i(triangle[e][2], triangle[e][0], triangle[e][1], node, x, y);
		}
	}

	return temp;
}

//���ģ�����
void error_max_2(int i, vector<double> cof_u, vector<vector<int>> triangle, vector<vector<double>> node, vector<int> node_boundary_marker)
{
	int m = 10, n = 10;
	double delta_x = 1.0 / m, delta_y = 1.0 / n;
	vector<vector<double>> error(m + 1, vector<double>(n + 1, 0));
	ofstream outfile;//�����ά��ֵ��
	string temp = "E:\\study_materials\\FEM\\Pg4\\matrix";
	string path = temp + to_string(i) + ".txt";
	cout << "���������������λ�ã�" << path << endl;
	outfile.open(path, ios::out);
	for (int i = 0; i <= m; i++)
		for (int j = 0; j <= n; j++)
			error[i][j] = fabs(u(i*delta_x, j*delta_y) - uh(i*delta_x, j*delta_y, cof_u, triangle, node, node_boundary_marker)),
			outfile << " " << i * delta_x << " " << j * delta_y << " " << error[i][j] << endl;
	outfile.close();
}

//���������
void matrix_output(vector<vector<double>> matrix)
{
	cout << endl << "�����������" << endl;
	for (int i = 0; i < matrix.size(); i++)
	{
		double temp = 0;
		cout << "{";
		for (int j = 0; j < matrix[0].size(); j++)
			cout << matrix[i][j] << ",", temp += matrix[i][j];
		cout << "},";
		temp -= matrix[i][i];
		cout << temp << endl;
	}
	cout << endl;
}

//���������
void vector_output(vector<double> vect)
{
	cout << endl << "�����������" << endl;
	for (int i = 0; i < vect.size(); i++)
		cout << vect[i] << ",";
	cout << endl;
}