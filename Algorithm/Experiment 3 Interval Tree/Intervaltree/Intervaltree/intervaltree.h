#pragma once

#include<iostream>
#include<vector>
#include<math.h>
#include <fstream>  //�ļ����⺯��

using namespace std;

//colorΪint�ͣ�0Ϊ��ɫ��1Ϊ��ɫ��
#define red 1
#define black 0
#define error -1

typedef struct TNode
{
	int color;
	int key;           //���������low
	int max;
	int low, high;
	struct TNode* left;
	struct TNode* right;
	struct TNode* parent;
}TNode;

typedef struct RBTree
{
	TNode* root;
	TNode* nil;
}RBTree;

//������Ϣ��ά��������maxֵ��
int attach(TNode* x)
{
	int temp = x->high;
	if (x->max != error && x->left->max > temp)
		temp = x->left->max;
	if (x->max != error && x->right->max > temp)
		temp = x->right->max;
	return temp;
}

//������
void LeftRotate(RBTree* T, TNode* x)
{
	TNode* y = x->right;
	x->right = y->left;
	x->max = attach(x);
	
	if (y->left != T->nil)
		y->left->parent = x;
	y->parent = x->parent;
	if (x->parent == T->nil)
		T->root = y;
	else if (x == x->parent->left)
		x->parent->left = y;
	else
		x->parent->right = y;
	y->parent->max = attach(y->parent);
	
	y->left = x;
	y->max = attach(y);

	x->parent = y;
}

//������
void RightRotate(RBTree* T, TNode* y)
{
	TNode* x = y->left;
	y->left = x->right;
	y->max = attach(y);

	if (x->right != T->nil)
		x->right->parent = y;
	x->parent = y->parent;
	if (y->parent == T->nil)
		T->root = x;
	else if (y == y->parent->right)
		y->parent->right = x;
	else
		y->parent->left = x;
	x->parent->max = attach(x->parent);

	x->right = y;
	x->max = attach(x);

	y->parent = x;
}

//�����������
void RBInsertFixup(RBTree* T, TNode* z)
{
	while (z->parent->color == red)  //��zΪ������z->parent=T->nil����ɫ�ڣ��������ѭ����
	{                                //��z->parentΪ�ڣ������������ͬ��������ѭ����
		if (z->parent == z->parent->parent->left)//case1,2,3
		{
			TNode* y = z->parent->parent->right; //y��z�����壻
			if (y->color == red) //case1
			{
				z->parent->color = black;
				y->color = black;
				z->parent->parent->color = red;
				z = z->parent->parent;
			}
			else //case2 or case3��yΪ��
			{
				if (z == z->parent->right) //case2
				{
					z = z->parent;
					LeftRotate(T, z);
				}//����Ϊcase3
				z->parent->color = black;
				z->parent->parent->color = red;
				RightRotate(T, z->parent->parent);
			}
		}
		else
		{
			TNode* y = z->parent->parent->left;//y��z�����壻
			if (y->color == red)//case4
			{
				z->parent->color = black;
				y->color = black;
				z->parent->parent->color = red;
				z = z->parent->parent;
			}
			else //case5 or case6��yΪ��
			{
				if (z == z->parent->left) //case5
				{
					z = z->parent;
					RightRotate(T, z);
				}//����Ϊcase6
				z->parent->color = black;
				z->parent->parent->color = red;
				LeftRotate(T, z->parent->parent);
			}
		}
	}
	T->root->color = black;
}

//���룻
void RBInsert(RBTree* T, int lowk, int highk)
{
	TNode* y = T->nil; //y��¼��ǰɨ��ڵ��˫�׽�㣻
	TNode* x = T->root;//�Ӹ���ʼɨ��
	TNode* z = (TNode*)malloc(sizeof(TNode));
	z->color = red;
	z->low = lowk;
	z->high = highk;
	z->key = lowk;
	z->max = highk;
	z->left = T->nil;
	z->right = T->nil;
	while (x != T->nil)
	{
		y = x;
		if (z->key < x->key)
			x = x->left;
		else
			x = x->right;
	}
	z->parent = y;     //y��z��˫��
	if (y == T->nil)   //z�����������z�Ǹ���
		T->root = z;
	else if (z->key < y->key)
		y->left = z;
	else
		y->right = z;
	RBInsertFixup(T, z);
}

//��ʼ���������
void InitRBTree(RBTree* T)
{
	//T = (RBTree *)malloc(sizeof(RBTree));
	//T->root= (TNode*)malloc(sizeof(TNode));
	T->nil = (TNode*)malloc(sizeof(TNode));
	T->nil->color = black;
	T->nil->key = error;
	T->nil->low = error;
	T->nil->high = error;
	T->nil->max = error;
	T->root = T->nil;
}

//�ж��Ƿ��ص������ص��򷵻�False��
bool overlap(TNode* x, int low, int high)
{
	if (x->low <= high && low <= x->high)
		return false;
	return true;
}

//�����ص����䣻
TNode* search(RBTree* T, int low, int high)
{
	TNode* x = T->root;
	while (x != T->nil && overlap(x, low, high))
	{
		if (x->left != T->nil && x->left->max >= low)
			x = x->left;
		else
			x = x->right;
	}
	return x;
}
