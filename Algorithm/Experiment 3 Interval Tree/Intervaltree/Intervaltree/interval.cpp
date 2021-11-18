#include"intervaltree.h"

using namespace std;

int main()
{
    vector<vector<int>> Re; //数组储存
    vector<int> temp(2, 0); //读入
    int n = 0;              //数组长度  
    RBTree* T;              //构造一棵树
    T = (RBTree*)malloc(sizeof(RBTree));

    ifstream infile;        //输入流
    //读入数据
    infile.open("E:\\study_materials\\Algorithm\\Intervaltree\\insert.txt", ios::in);
    if (!infile.is_open())
        cout << "Open file failure" << endl;
    infile >> n;
    for (int i = 0; i < n; i++)
    {
        infile >> temp[0];
        infile >> temp[1];
        Re.push_back(temp);
    }
    infile.close();
    //构建区间树；
    InitRBTree(T);
    for (int i = 0; i < n; i++)
        RBInsert(T, Re[i][0], Re[i][1]);

    int request = 1;
    int low = 0, high = 0;
    TNode* interval=(TNode*)malloc(sizeof(TNode));
    while (request != 0)
    {
        cout << endl << "请输入待查找的区间（先输入低点，再输入高点）:" << endl;
        cin >> low >> high;
        interval = search(T, low, high);
        if (interval == T->nil)
        {
            cout << "查找区间为：[" << low << "," << high << "]；   ";
            cout << "无重叠区间" << endl << endl;
        }
        else
        {
            cout << "查找区间为：[" << low << "," << high << "]；   ";
            cout << "重叠区间为：[" << interval->low << "," << interval->high << "]" << endl << endl;
        }
        cout << "是否继续查找，若否请输入0，若继续查找请输入1：" << endl;
        cin >> request;
    }
}