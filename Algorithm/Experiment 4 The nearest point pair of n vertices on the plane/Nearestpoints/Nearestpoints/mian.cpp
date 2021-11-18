#include <iostream>
#include <cmath>

using namespace std;

//点结构；
struct point
{
	double x, y;
	char  num;
};
//最近点对结构；
struct closestpoint
{
	point left, right;
	double distance;
};
//平面点对；
point points[6] = { {1.1,2.0,'A'},{3.1,1.0,'B'},
	{4.4,5.3,'C'},{2.6,3.2,'D'},
	{2.6,3.4,'E'},{5.3,5.2,'F'} };
point midpoint[6];

//依据横坐标排序，用于初始时排序；
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

//依据纵坐标排序，用于递归中排序；
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

//距离公式；
double dist(int i, int j)
{
	return sqrt((points[i].x - points[j].x) * (points[i].x - points[j].x) + (points[i].y - points[j].y) * (points[i].y - points[j].y));
}

//求最近点对；
closestpoint nearestpoints(int low, int high)
{
	closestpoint temp = { points[low],points[high],dist(low,high) };
	switch (high - low + 1)
	{
	case 1: return temp; //一个点默认返回最近点对均为自己。
	case 2: return temp; //两个点就返回这两个点即可。
	case 3: //三个点暴力求解。
	{
		if (dist(low, low + 1) < temp.distance)
			temp = { points[low],points[low + 1],dist(low,low + 1) };
		if (dist(low + 1, high) < temp.distance)
			temp = { points[low + 1],points[high],dist(low + 1,high) };
		return temp;
	}
	default: //>3的情况
	{
		closestpoint left = nearestpoints(low, floor((low + high) / 2));
		closestpoint right = nearestpoints(floor((low + high) / 2) + 1, high);
		temp = (left.distance <= right.distance) ? left : right;//取得左右最小中更小的点对；
		
		//point midpoint[size(points)];
		int length = 0;
		int m = floor((low + high) / 2);
		for (int i = 0; i < size(points); i++)
			if (abs(points[i].x - points[m].x) <= temp.distance)
			{
				midpoint[length] = points[i];
				length++;
			}//所以处在中间待判断的共length个；

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

//主函数；
int main()
{
	cout << "平面点对为：" << endl;
	for (int i = 0; i < size(points); i++)
		cout << points[i].num << " (" << points[i].x << "," << points[i].y << ")" << endl;

	cout << endl;

	//排序；
	insertsort(0);

	closestpoint min;
	min = nearestpoints(0, size(points) - 1);

	cout << "最近点对为：" << min.left.num << " " << min.right.num << endl;
	cout << "最近距离为：" << min.distance << endl;
}
