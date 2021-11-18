#pragma once

#include<iostream>
#include<vector>
#include<math.h>
#include <fstream>  //文件流库函数

using namespace std;

//color为int型，0为黑色，1为红色；
#define red 1
#define black 0
#define error -1

typedef struct TNode
{
	int color;
	int key;           //储存区间的low
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

//附加信息的维护，返回max值；
int attach(TNode* x)
{
	int temp = x->high;
	if (x->max != error && x->left->max > temp)
		temp = x->left->max;
	if (x->max != error && x->right->max > temp)
		temp = x->right->max;
	return temp;
}

//左旋；
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

//右旋；
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

//插入的修正；
void RBInsertFixup(RBTree* T, TNode* z)
{
	while (z->parent->color == red)  //若z为根，则z->parent=T->nil，颜色黑，不进入此循环；
	{                                //若z->parent为黑，则无需调整，同样不进入循环；
		if (z->parent == z->parent->parent->left)//case1,2,3
		{
			TNode* y = z->parent->parent->right; //y是z的叔叔；
			if (y->color == red) //case1
			{
				z->parent->color = black;
				y->color = black;
				z->parent->parent->color = red;
				z = z->parent->parent;
			}
			else //case2 or case3，y为黑
			{
				if (z == z->parent->right) //case2
				{
					z = z->parent;
					LeftRotate(T, z);
				}//以下为case3
				z->parent->color = black;
				z->parent->parent->color = red;
				RightRotate(T, z->parent->parent);
			}
		}
		else
		{
			TNode* y = z->parent->parent->left;//y是z的叔叔；
			if (y->color == red)//case4
			{
				z->parent->color = black;
				y->color = black;
				z->parent->parent->color = red;
				z = z->parent->parent;
			}
			else //case5 or case6，y为黑
			{
				if (z == z->parent->left) //case5
				{
					z = z->parent;
					RightRotate(T, z);
				}//以下为case6
				z->parent->color = black;
				z->parent->parent->color = red;
				LeftRotate(T, z->parent->parent);
			}
		}
	}
	T->root->color = black;
}

//插入；
void RBInsert(RBTree* T, int lowk, int highk)
{
	TNode* y = T->nil; //y记录当前扫描节点的双亲结点；
	TNode* x = T->root;//从根开始扫描
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
	z->parent = y;     //y是z的双亲
	if (y == T->nil)   //z插入空树，故z是根；
		T->root = z;
	else if (z->key < y->key)
		y->left = z;
	else
		y->right = z;
	RBInsertFixup(T, z);
}

//初始化红黑树；
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

//判定是否重叠，若重叠则返回False；
bool overlap(TNode* x, int low, int high)
{
	if (x->low <= high && low <= x->high)
		return false;
	return true;
}

//查找重叠区间；
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
