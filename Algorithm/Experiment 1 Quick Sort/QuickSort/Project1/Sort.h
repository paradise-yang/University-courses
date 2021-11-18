#pragma once

#include<iostream>
#include<vector>
#include<math.h>

using namespace std;

//���ţ���׼ΪA[r];
int partition_last(vector<int>& A, int p, int r)
{
	int x = A[r];
	int i = p - 1;
	int temp = 0;
	for (int j = p; j <= r - 1; j++)
	{
		if (A[j] <= x)
		{
			i = i + 1;
			temp = A[i];
			A[i] = A[j];
			A[j] = temp;
		}
	}
	temp = A[i + 1];
	A[i + 1] = A[r];
	A[r] = temp;
	return i + 1;
}

//��׼����;
void quicksort(vector<int>& A, int p, int r)
{
	if (p < r)
	{
		int q = partition_last(A, p, r);
		quicksort(A, p, q - 1);
		quicksort(A, q + 1, r);
	}
}

//���ţ���׼���;
int randomized_partition(vector<int>& A, int p, int r)
{
	int i = (rand() % (r - p + 1)) + p; //����[p,r]����������
	int temp = A[r];
	A[r] = A[i];
	A[i] = temp;
	return partition_last(A,p,r);
}

//������ţ�
void quicksort_random(vector<int>& A, int p, int r)
{
	if (p < r)
	{
		int q = randomized_partition(A, p, r);
		quicksort(A, p, q - 1);
		quicksort(A, q + 1, r);
	}
}

//����ȡ�л�׼��
int partition_mid(vector<int>& A, int p, int r)
{
	int a = A[p], b = A[floor((p + r) / 2)], c = A[r];
	int x = 0;

	if ((a - b) * (a - c) <= 0)
		x = a;
	if ((b - a) * (b - c) <= 0)
		x = b;
	if ((c - a) * (c - b) <= 0)
		x = c;

	int i = p - 1;
	int temp = 0;
	for (int j = p; j <= r - 1; j++)
	{
		if (A[j] <= x)
		{
			i = i + 1;
			temp = A[i];
			A[i] = A[j];
			A[j] = temp;
		}
	}
	temp = A[i + 1];
	A[i + 1] = A[r];
	A[r] = temp;
	return i + 1;
}

//����ȡ�п��ţ�
void quicksort_mid(vector<int>& A, int p, int r)
{
	if (p < r)
	{
		int q = partition_mid(A, p, r);
		quicksort(A, p, q - 1);
		quicksort(A, q + 1, r);
	}
}

//��������
void insertsort(vector<int>& A)
{
	int key = 0, i = 0;
	for (int j = 1; j < A.size(); j++)
	{
		key = A[j];
		i = j - 1;
		while (i >= 0 && A[i] > key)
		{
			A[i + 1] = A[i];
			i = i - 1;
		}
		A[i + 1] = key;
	}
}

//���䳤��<=kʱ�������ţ�ֱ�ӷ��أ�
void quicksort_k(vector<int>& A, int p, int r, int k)
{
	if (p<r && r - p >= k)
	{
		int q = partition_last(A, p, r);
		quicksort_k(A, p, q - 1, k);
		quicksort_k(A, q + 1, r, k);
	}
}

//��������Ŀ��ţ�
void quicksort_insert(vector<int>& A, int p, int r, int k)
{
	quicksort_k(A, p, r, k);
	//for (int i = 0; i < A.size(); i++)
	//{
	//	cout << "  " << A[i];
	//}
	insertsort(A);
}
