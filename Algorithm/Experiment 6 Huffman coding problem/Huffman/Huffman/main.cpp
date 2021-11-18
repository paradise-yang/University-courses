#include<iostream>
#include<vector>
#include<math.h>
#include<string.h>
#include<stdio.h>
#include<string>
#include<map>
#include <queue>
#include<algorithm>
#include <fstream>  //文件流库函数

#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable:4996);

using namespace std;
#define MAX_LINE 1024
#define LEN 512
struct huffman_node {
    char c;
    int weight;
    char huffman_code[LEN];
    huffman_node* left;
    huffman_node* right;
};

//读取文件并解析词频
int read_file(FILE* fn, map<char, int>& word, map<char, string>& code)
{
    if (fn == NULL) return 0;
    char line[MAX_LINE];
    while (fgets(line, MAX_LINE, fn)) 
    {//解析，统计词频
        char* p = line;
        while (*p != '\0' && *p != '\n') {
            map<char, int>::iterator it = word.find(*p);
            if (it == word.end()) // 不存在，插入
                code.insert(make_pair(*p, "\0")), word.insert(make_pair(*p, 1));
            else
                it->second++;
            p++;
        }
    }
    return 0;
}

//按weight升序排列
bool sort_by_weight(huffman_node* a, huffman_node* b) { return (a->weight < b->weight); }

//构造Huffman树
int huffman_tree_create(huffman_node*& root, map<char, int>& word) {
    char line[MAX_LINE];
    vector<huffman_node*> huffman_tree_node;
    map<char, int>::iterator it_t;
    for (it_t = word.begin(); it_t != word.end(); it_t++) // 为每一个节点申请空间
    {
        huffman_node* node = (huffman_node*)malloc(sizeof(huffman_node));
        node->c = it_t->first;
        node->weight = it_t->second;
        node->huffman_code[0] = '\0';
        node->left = NULL;
        node->right = NULL;
        huffman_tree_node.push_back(node);
    }
    while (huffman_tree_node.size() > 0) //从叶节点开始构建Huffman树
    {
        sort(huffman_tree_node.begin(), huffman_tree_node.end(), sort_by_weight);// 按照weight升序排序
        /*
        for (int i = 0; i < huffman_tree_node.size(); i++)
            cout << huffman_tree_node[i]->c << " " << huffman_tree_node[i]->weight << endl;
        cout << endl;
        */
        if (huffman_tree_node.size() == 1) // 只有一个根结点
        {
            root = huffman_tree_node[0];
            huffman_tree_node.erase(huffman_tree_node.begin());
        }
        else 
        {   // 取出前两个
            huffman_node* node_1 = huffman_tree_node[0];
            huffman_node* node_2 = huffman_tree_node[1];
            // 删除
            huffman_tree_node.erase(huffman_tree_node.begin());
            huffman_tree_node.erase(huffman_tree_node.begin());
            // 生成新的节点
            huffman_node* node = (huffman_node*)malloc(sizeof(huffman_node));
            node->c = '.', node->left = NULL, node->right = NULL, node->huffman_code[0] = '\0';
            node->weight = node_1->weight + node_2->weight;
            (node_1->weight < node_2->weight) ? (node->left = node_1, node->right = node_2) : (node->left = node_2, node->right = node_1);
            huffman_tree_node.push_back(node);
        }
    }
    return 0;
}
/*
//先序遍历，只输出带字母的节点
void print_huffman_pre(huffman_node* node) {
    if (node != NULL) {
        if (node->c != '.')
            cout << node->c << " " << node->weight << endl;
        print_huffman_pre(node->left);
        print_huffman_pre(node->right);
    }
}
//中序遍历，只输出带字母的节点
void print_huffman_in(huffman_node* node) {
    if (node != NULL) {
        print_huffman_in(node->left);
        if (node->c != '.')
            cout << node->c << " " << node->weight << endl;
        print_huffman_in(node->right);
    }
}
*/
//实现Huffman编码
int get_huffman_code(huffman_node*& node) 
{
    if (node == NULL) return 1;
    // 利用层次遍历，构造每一个节点
    huffman_node* p = node;
    queue<huffman_node*> q;
    q.push(p);
    while (q.size() > 0) {
        p = q.front();
        q.pop();
        if (p->left != NULL) 
        {
            q.push(p->left);
            strcpy((p->left)->huffman_code, p->huffman_code);//, size(p->huffman_code)
            char* ptr = (p->left)->huffman_code;
            while (*ptr != '\0') ptr++;
            *ptr = '0';
            ptr++;
            *ptr = '\0';
        }
        if (p->right != NULL) 
        {
            q.push(p->right);
            strcpy((p->right)->huffman_code, p->huffman_code);
            char* ptr = (p->right)->huffman_code;
            while (*ptr != '\0') ptr++;
            *ptr = '1';
            ptr++;
            *ptr = '\0';
        }
    }
    return 0;
}

//打印叶子结点信息，并计算对应长度
void print_leaf(huffman_node* node, map<char, int>& word, int& huffman_size, map<char, string>& code)
{
    if (node != NULL) 
    {
        string temp;
        print_leaf(node->left, word, huffman_size, code);
        if (node->left == NULL && node->right == NULL)
        {
            temp = node->huffman_code;
            map<char, int>::iterator it = word.find(node->c);
            map<char, string>::iterator its = code.find(node->c);
            its->second = temp;
            huffman_size += size(temp) * it->second;
            //cout << node->c << " " << node->huffman_code << " " << size(temp) << endl;
        }
        print_leaf(node->right, word, huffman_size, code);
    }
}

//计算压缩率
void compression_ratio(map<char, int> word, int huffman_size)
{
    //计算定长编码的长度
    int k = 0, n = word.size() - 1;
    while (n)
    {
        ++k;
        n >>= 1;
    }
    //计算字符数
    n = 0;
    map<char, int>::iterator m1_Iter;
    for (m1_Iter = word.begin(); m1_Iter != word.end(); m1_Iter++)
        n += m1_Iter->second;
    cout << "Huffman编码所需字符长度为：" << huffman_size << endl;
    cout << "定长编码所需字符长度为：" << k * n << endl;
    cout << "压缩率为：" << 1.0 * huffman_size / (k * n) * 100 << "%" << endl;

}

//释放Huffman树
void destory_huffman_tree(huffman_node* node) 
{
    if (node != NULL) 
    {
        destory_huffman_tree(node->left);
        destory_huffman_tree(node->right);
        free(node);
        node = NULL;
    }
}

//输出Huffman编码到指定文件
int output(FILE* fn, ofstream& outfile, map<char, string>& code)
{
    if (fn == NULL) return 0;
    char line[MAX_LINE];
    while (fgets(line, MAX_LINE, fn))
    {
        char* p = line;
        while (*p != '\0' && *p != '\n') {
            map<char, string>::iterator it = code.find(*p);
            outfile << it->second;
            p++;
        }
    }
    return 0;
}

//读取文件路径
#define F_PATH "E:\\study_materials\\Algorithm\\Huffman\\data.txt"

int main() {
    // 读文件
    FILE* fn = fopen(F_PATH, "r");
    huffman_node* root = NULL;
    map<char, int> word;
    map<char, string>code;
    read_file(fn, word, code);
    /*
    map<char, int>::iterator m1_Iter;
    for (m1_Iter = word.begin(); m1_Iter != word.end(); m1_Iter++)
        cout << m1_Iter->first << " " << m1_Iter->second << endl;
    */
    huffman_tree_create(root, word);
    fclose(fn);
    /*
    cout << "in-order:" << endl;
    print_huffman_in(root);
    */
    get_huffman_code(root);
    cout << "Huffman编码为：" << endl;
    int huffman_size = 0;
    print_leaf(root, word, huffman_size, code);
    map<char, string>::iterator m1_Iter;
    for (m1_Iter = code.begin(); m1_Iter != code.end(); m1_Iter++)
        cout << m1_Iter->first << " " << m1_Iter->second << endl;
    cout << endl;
    compression_ratio(word, huffman_size);

    fn = fopen(F_PATH, "r");
    ofstream outfile;       //输出流
    outfile.open("E:\\study_materials\\Algorithm\\Huffman\\encode.txt", ios::out);
    output(fn, outfile, code);
    outfile.close();
    fclose(fn);
    destory_huffman_tree(root);
}
