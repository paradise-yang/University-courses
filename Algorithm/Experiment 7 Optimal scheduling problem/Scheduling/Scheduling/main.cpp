#include<iostream>

using namespace std;

#define NUM_TASK 10
#define NUM_MACHSHINE 4

int x[NUM_TASK] = { 0 };         //��¼���䣬��x[task]��ʾ����task���������x[task]
int best_x[NUM_TASK] = { 0 };    //�洢���ŷ��䷽��
double min_time = 1024;          //ִ������������Сʱ��
double time_task[NUM_TASK] = { 1.0,14.3,10.3,2.4,5.3,3.1,4.9,2.1,8.4,15.9 };//ÿ����������ʱ��
double time_machine_end[NUM_MACHSHINE] = { 0.0 };                            //ÿ���������н���ʱ��

//��ȡ��ǰ�ѷ�����������ʱ��
double getmaxtime(double time_mac[]) 
{
	double max_time = time_mac[0];
	for (int i = 1; i < NUM_MACHSHINE; i++) 
		if (time_mac[i] > max_time) 
			max_time = time_mac[i];
	return max_time;
}

//���ݷ����
void BackTrack(int task) 
{
	if (task >= NUM_TASK)
	{
		double current_time = getmaxtime(time_machine_end);//��ǰ�ѷ�����������ʱ��
		if (current_time < min_time) 
		{
			min_time = current_time;
			for (int i = 0; i < NUM_TASK; i++)
				best_x[i] = x[i];
		}
	}
	else 
	{
		for (int i = 0; i < NUM_MACHSHINE; i++) 
		{
			x[task] = i;
			time_machine_end[i] += time_task[task];
			if (time_machine_end[i] < min_time)
				BackTrack(task + 1);
			time_machine_end[i] -= time_task[task];
		}
	}
}

//������䷽��
void output_assign(int best_x[])
{
	for (int i = 0; i < NUM_TASK; i++)
		cout << "����" << i + 1 << "���������" << best_x[i] + 1 << endl;
}

int main() {
	BackTrack(0);

	cout << "��������ִ��ʱ������Ϊ��";
	for (int i = 0; i < NUM_TASK; i++)
		cout << time_task[i] << " ";

	cout << endl << endl;
	cout << "���������������Ҫ����Сʱ��Ϊ��" << min_time;

	cout << endl << endl;
	output_assign(best_x);

	return 0;
}
