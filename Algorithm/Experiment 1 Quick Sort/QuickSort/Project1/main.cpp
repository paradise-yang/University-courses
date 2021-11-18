#include<iostream>
#include <vector>
#include <fstream>  //文件流库函数
#include <chrono>
//#include <string>
#include "Sort.h"

using namespace std;

int main()
{
    vector<int> Re;         //数组储存
    vector<int> temp(1,0);  //读入
    int n = 0;              //数组长度   
    ifstream infile;        //输入流
    ofstream outfile;       //输出流
    
    infile.open("E:\\study_materials\\Algorithm\\QuickSort\\data.txt", ios::in);
    if (!infile.is_open())
        cout << "Open file failure" << endl;
    infile >> n;
    for (int i = 0; i < n; i++)
    {
        infile >> temp[0];
        Re.push_back(temp[0]);
    }
        infile.close();   //关闭data文件

    vector<int> Backup = Re;
    cout << "固定基准快排所需时间:";
    auto start = std::chrono::steady_clock::now();
    quicksort(Re, 0, n - 1);
    auto end = chrono::steady_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << elapsed_seconds.count() << "s" << endl;
    
    outfile.open("E:\\study_materials\\Algorithm\\QuickSort\\sorted.txt", ios::app);
    if (!outfile.is_open())
        cout << "Open file failure" << endl;
    for (int i = 0; i < n; i++)
        outfile << Re[i] << " ";
    outfile.close();

    Re = Backup;
    cout << "随机基准快排所需时间:";
    start = std::chrono::steady_clock::now();
    quicksort_random(Re, 0, n - 1);
    end = chrono::steady_clock::now();
    elapsed_seconds = end - start;
    cout << elapsed_seconds.count() << "s" << endl;


    Re = Backup;
    cout << "三数取中快排所需时间:";
    start = std::chrono::steady_clock::now();
    quicksort_mid(Re, 0, n - 1);
    end = chrono::steady_clock::now();
    elapsed_seconds = end - start;
    cout << elapsed_seconds.count() << "s" << endl;
    

    Re = Backup;
    cout << "快排结合插入所需时间:";
    start = std::chrono::steady_clock::now();
    quicksort_insert(Re, 0, n - 1, 3);
    end = chrono::steady_clock::now();
    elapsed_seconds = end - start;
    cout << elapsed_seconds.count() << "s" << endl;
    

    Re = Backup;
    cout << "对比插入排序所需时间:";
    start = std::chrono::steady_clock::now();
    insertsort(Re);
    end = chrono::steady_clock::now();
    elapsed_seconds = end - start;
    cout << elapsed_seconds.count() << "s" << endl;
    

    return 0;
    while (1);

}
