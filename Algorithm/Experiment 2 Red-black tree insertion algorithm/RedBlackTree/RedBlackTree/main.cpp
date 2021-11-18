#include"redblack.h"

using namespace std;

int main()
{
    vector<int> Re;         //���鴢��
    vector<int> temp(1, 0); //����
    int n = 0;              //���鳤��  
    RBTree *T;              //����һ����
    T = (RBTree*)malloc(sizeof(RBTree));
    
    ifstream infile;        //������
    ofstream outfile;       //�����

    infile.open("E:\\study_materials\\Algorithm\\RedBlackTree\\insert.txt", ios::in);
    if (!infile.is_open())
        cout << "Open file failure" << endl;
    infile >> n;
    for (int i = 0; i < n; i++)
    {
        infile >> temp[0];
        Re.push_back(temp[0]);
    }
    infile.close();

    InitRBTree(T);
    for (int i = 0; i < n; i++)
        RBInsert(T, Re[i]);

    outfile.open("E:\\study_materials\\Algorithm\\RedBlackTree\\LNR.txt", ios::trunc);
    if (!outfile.is_open())
        cout << "Open file failure" << endl;
    //cout << "������������" << endl;
    LNR(T->root, outfile);
    outfile.close();

    outfile.open("E:\\study_materials\\Algorithm\\RedBlackTree\\NLR.txt", ios::trunc);
    if (!outfile.is_open())
        cout << "Open file failure" << endl;
    //cout << "������������" << endl;
    NLR(T->root, outfile);
    outfile.close();
}