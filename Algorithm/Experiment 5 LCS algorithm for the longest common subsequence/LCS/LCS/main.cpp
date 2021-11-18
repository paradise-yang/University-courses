#include<iostream>
#include<vector>
#include<math.h>
#include<string.h>
#include<stdio.h>
#include<string>

using namespace std;

#define left 1
#define up 2
#define upleft 3

//求最大公共子序列；时间复杂度为O(mn)，空间复杂度为O(mn)
void LCS_mn(string x, string y, vector<vector<int>>& b, vector<vector<int>>& c)
{
	for (int i = 1; i <= size(x); i++)
	{
		for (int j = 1; j <= size(y); j++)
		{
			if (x[i - 1] == y[j - 1])
				c[i][j] = c[i - 1][j - 1] + 1, b[i][j] = upleft;
			else if (c[i - 1][j] >= c[i][j - 1])
				c[i][j] = c[i - 1][j], b[i][j] = up;
			else
				c[i][j] = c[i][j - 1], b[i][j] = left;
		}
	}
}

//输出最大公共子序列；
char printLCS(vector<vector<int>> b, string x, int i, int j)
{
	if (i == 0 || j == 0)
		return ' ';
	if (b[i][j] == upleft)
		printLCS(b, x, i - 1, j - 1), cout << x[i - 1];
	else if (b[i][j] == up)
		printLCS(b, x, i - 1, j);
	else
		printLCS(b, x, i, j - 1);
	
}

//求最大公共子序列长度；时间复杂度为O(mn)，空间复杂度为O(2min{m,n})
int LCS_2n(string x, string y)
{
	if (size(x) < size(y))//取最小的字符串放在y里面；
	{
		string temp = y;
		y = x;
		x = temp;
	}
	vector<vector<int>> c(2, vector<int>(size(y) + 1, 0));
	for (int i = 1; i <= size(x); i++)
	{
		for (int j = 1; j <= size(y); j++)
		{
			if (x[i - 1] == y[j - 1])
				c[i % 2][j] = c[(i + 1) % 2][j - 1] + 1;
			else if (c[(i + 1) % 2][j] >= c[i % 2][j - 1])
				c[i % 2][j] = c[(i + 1) % 2][j];
			else
				c[i % 2][j] = c[i % 2][j - 1];
		}
	}
	return c[size(x) % 2][size(y)];
}

//求最大公共子序列长度；时间复杂度为O(mn)，空间复杂度为O(min{m,n})
int LCS_n(string x, string y)
{
	if (size(x) < size(y))//取最小的字符串放在y里面；
	{
		string temp = y;
		y = x;
		x = temp;
	}
	vector<int> c(size(y) + 1, 0);
	for (int i = 1; i < size(x) + 1; i++)
	{
		int temp = c[0];
		for (int j = 1; j < size(y) + 1; j++)
		{
			if (x[i - 1] == y[j - 1])
			{
				int tempp = c[j];
				c[j] = temp + 1;
				temp = tempp;
			}
			else if (c[j] >= c[j - 1])
				temp = c[j];
			else
				temp = c[j], c[j] = c[j - 1];
		}
	}
	return c[size(y)];
}

int main()
{
	int require = 1;
	while (require)
	{
		cout << "请输入第一个字符串：";
		string text1;
		getline(cin, text1);
		cout << "请输入第二个字符串：";
		string text2;
		getline(cin, text2);

		//空间复杂度为O(mn);
		vector<vector<int>> b(size(text1) + 1, vector<int>(size(text2) + 1, 0));
		vector<vector<int>> c(size(text1) + 1, vector<int>(size(text2) + 1, 0));
		LCS_mn(text1, text2, b, c);
		cout << endl << "空间复杂度为O(mn)的算法，返回LCS的长度与序列：" << endl;
		if (c[size(text1)][size(text2)] != 0)
		{
			cout << "LCS：";
			printLCS(b, text1, size(text1), size(text2));
			cout << endl << "长度：" << c[size(text1)][size(text2)] << endl;
		}
		else
			cout << "长度：" << 0 << endl;

		//空间复杂度为O(2*min{m,n});
		cout << endl << "空间复杂度为O(2*min{m,n})的算法，返回LCS的长度与序列：";
		cout << endl << "长度：" << LCS_2n(text1, text2) << endl;

		//空间复杂度为O(min{m,n});
		cout << endl << "空间复杂度为O(min{m,n})的算法，返回LCS的长度与序列：";
		cout << endl << "长度：" << LCS_n(text1, text2) << endl;

		cout << endl << endl;
		cout << endl << "是否继续查找LCS，若是输入任何非0数，若否输入0结束：";
		cin >> require;
		cout << endl;
		getchar();
	}
}