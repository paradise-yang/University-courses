#include"intervaltree.h"

using namespace std;

int main()
{
    vector<vector<int>> Re; //���鴢��
    vector<int> temp(2, 0); //����
    int n = 0;              //���鳤��  
    RBTree* T;              //����һ����
    T = (RBTree*)malloc(sizeof(RBTree));

    ifstream infile;        //������
    //��������
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
    //������������
    InitRBTree(T);
    for (int i = 0; i < n; i++)
        RBInsert(T, Re[i][0], Re[i][1]);

    int request = 1;
    int low = 0, high = 0;
    TNode* interval=(TNode*)malloc(sizeof(TNode));
    while (request != 0)
    {
        cout << endl << "����������ҵ����䣨������͵㣬������ߵ㣩:" << endl;
        cin >> low >> high;
        interval = search(T, low, high);
        if (interval == T->nil)
        {
            cout << "��������Ϊ��[" << low << "," << high << "]��   ";
            cout << "���ص�����" << endl << endl;
        }
        else
        {
            cout << "��������Ϊ��[" << low << "," << high << "]��   ";
            cout << "�ص�����Ϊ��[" << interval->low << "," << interval->high << "]" << endl << endl;
        }
        cout << "�Ƿ�������ң�����������0������������������1��" << endl;
        cin >> request;
    }
}