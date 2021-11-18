#include <iostream>
#include <cmath>

using namespace std;

//��ṹ��
struct point
{
	double x, y;
	char  num;
};
//�����Խṹ��
struct closestpoint
{
	point left, right;
	double distance;
};
//ƽ���ԣ�
point points[6] = { {1.1,2.0,'A'},{3.1,1.0,'B'},
	{4.4,5.3,'C'},{2.6,3.2,'D'},
	{2.6,3.4,'E'},{5.3,5.2,'F'} };
point midpoint[6];

//���ݺ������������ڳ�ʼʱ����
void insertsort(int sort)
{
	point key = { 0,0,'o' };
	int i = 0;
	for (int j = 1; j < size(points); j++)
	{
		key = points[j];
		i = j - 1;
		while (i >= 0 && points[i].x > key.x)
		{
			points[i + 1] = points[i];
			i = i - 1;
		}
		points[i + 1] = key;
	}
}

//�����������������ڵݹ�������
void sorty(int end)
{
	point key = { 0,0,'o' };
	int i = 0;
	for (int j = 1; j < end; j++)
	{
		key = midpoint[j];
		i = j - 1;
		while (i >= 0 && midpoint[i].y > key.y)
		{
			midpoint[i + 1] = midpoint[i];
			i = i - 1;
		}
		midpoint[i + 1] = key;
	}
}

//���빫ʽ��
double dist(int i, int j)
{
	return sqrt((points[i].x - points[j].x) * (points[i].x - points[j].x) + (points[i].y - points[j].y) * (points[i].y - points[j].y));
}

//�������ԣ�
closestpoint nearestpoints(int low, int high)
{
	closestpoint temp = { points[low],points[high],dist(low,high) };
	switch (high - low + 1)
	{
	case 1: return temp; //һ����Ĭ�Ϸ��������Ծ�Ϊ�Լ���
	case 2: return temp; //������ͷ����������㼴�ɡ�
	case 3: //�����㱩����⡣
	{
		if (dist(low, low + 1) < temp.distance)
			temp = { points[low],points[low + 1],dist(low,low + 1) };
		if (dist(low + 1, high) < temp.distance)
			temp = { points[low + 1],points[high],dist(low + 1,high) };
		return temp;
	}
	default: //>3�����
	{
		closestpoint left = nearestpoints(low, floor((low + high) / 2));
		closestpoint right = nearestpoints(floor((low + high) / 2) + 1, high);
		temp = (left.distance <= right.distance) ? left : right;//ȡ��������С�и�С�ĵ�ԣ�
		
		//point midpoint[size(points)];
		int length = 0;
		int m = floor((low + high) / 2);
		for (int i = 0; i < size(points); i++)
			if (abs(points[i].x - points[m].x) <= temp.distance)
			{
				midpoint[length] = points[i];
				length++;
			}//���Դ����м���жϵĹ�length����

		sorty(length);

		for(int i=0;i<length;i++)
			for (int j = i + 1; j < ((i + 6 < length) ? (i + 6) : length); j++)
			{
				if (dist(i, j) >= temp.distance)
					break;
				else
					temp = { points[i],points[j],dist(i,j) };
			}

		return temp;
	}
	}
}

//��������
int main()
{
	cout << "ƽ����Ϊ��" << endl;
	for (int i = 0; i < size(points); i++)
		cout << points[i].num << " (" << points[i].x << "," << points[i].y << ")" << endl;

	cout << endl;

	//����
	insertsort(0);

	closestpoint min;
	min = nearestpoints(0, size(points) - 1);

	cout << "������Ϊ��" << min.left.num << " " << min.right.num << endl;
	cout << "�������Ϊ��" << min.distance << endl;
}
