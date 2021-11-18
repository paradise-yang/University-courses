#include<iostream>
#include <vector>
#include <fstream>  //�ļ����⺯��
#include <chrono>
//#include <string>
#include "Sort.h"

using namespace std;

int main()
{
    vector<int> Re;         //���鴢��
    vector<int> temp(1,0);  //����
    int n = 0;              //���鳤��   
    ifstream infile;        //������
    ofstream outfile;       //�����
    
    infile.open("E:\\study_materials\\Algorithm\\QuickSort\\data.txt", ios::in);
    if (!infile.is_open())
        cout << "Open file failure" << endl;
    infile >> n;
    for (int i = 0; i < n; i++)
    {
        infile >> temp[0];
        Re.push_back(temp[0]);
    }
        infile.close();   //�ر�data�ļ�

    vector<int> Backup = Re;
    cout << "�̶���׼��������ʱ��:";
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
    cout << "�����׼��������ʱ��:";
    start = std::chrono::steady_clock::now();
    quicksort_random(Re, 0, n - 1);
    end = chrono::steady_clock::now();
    elapsed_seconds = end - start;
    cout << elapsed_seconds.count() << "s" << endl;


    Re = Backup;
    cout << "����ȡ�п�������ʱ��:";
    start = std::chrono::steady_clock::now();
    quicksort_mid(Re, 0, n - 1);
    end = chrono::steady_clock::now();
    elapsed_seconds = end - start;
    cout << elapsed_seconds.count() << "s" << endl;
    

    Re = Backup;
    cout << "���Ž�ϲ�������ʱ��:";
    start = std::chrono::steady_clock::now();
    quicksort_insert(Re, 0, n - 1, 3);
    end = chrono::steady_clock::now();
    elapsed_seconds = end - start;
    cout << elapsed_seconds.count() << "s" << endl;
    

    Re = Backup;
    cout << "�ԱȲ�����������ʱ��:";
    start = std::chrono::steady_clock::now();
    insertsort(Re);
    end = chrono::steady_clock::now();
    elapsed_seconds = end - start;
    cout << elapsed_seconds.count() << "s" << endl;
    

    return 0;
    while (1);

}
