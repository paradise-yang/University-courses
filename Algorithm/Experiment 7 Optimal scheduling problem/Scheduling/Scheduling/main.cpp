#include<iostream>

using namespace std;

#define NUM_TASK 10
#define NUM_MACHSHINE 4

int x[NUM_TASK] = { 0 };         //记录分配，即x[task]表示任务task分配给机器x[task]
int best_x[NUM_TASK] = { 0 };    //存储最优分配方案
double min_time = 1024;          //执行任务所需最小时间
double time_task[NUM_TASK] = { 1.0,14.3,10.3,2.4,5.3,3.1,4.9,2.1,8.4,15.9 };//每个任务所需时间
double time_machine_end[NUM_MACHSHINE] = { 0.0 };                            //每个机器运行结束时间

//获取当前已分配任务的完成时间
double getmaxtime(double time_mac[]) 
{
	double max_time = time_mac[0];
	for (int i = 1; i < NUM_MACHSHINE; i++) 
		if (time_mac[i] > max_time) 
			max_time = time_mac[i];
	return max_time;
}

//回溯法求解
void BackTrack(int task) 
{
	if (task >= NUM_TASK)
	{
		double current_time = getmaxtime(time_machine_end);//当前已分配任务的完成时间
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

//输出分配方案
void output_assign(int best_x[])
{
	for (int i = 0; i < NUM_TASK; i++)
		cout << "任务" << i + 1 << "分配给机器" << best_x[i] + 1 << endl;
}

int main() {
	BackTrack(0);

	cout << "各个任务执行时间依次为：";
	for (int i = 0; i < NUM_TASK; i++)
		cout << time_task[i] << " ";

	cout << endl << endl;
	cout << "所有任务完成所需要的最小时间为：" << min_time;

	cout << endl << endl;
	output_assign(best_x);

	return 0;
}
