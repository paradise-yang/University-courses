#include <iostream>
#include <queue>

using namespace std;

#define vexnumber  8
char symbol[vexnumber] = { 'A','B','C','D','E','F','G','H' };
int first[] = { 0,0,2,2,3,5 };
int second[] = { 1,2,3,4,4,6 };

typedef int infoType;
typedef char vertexType;
typedef int Status;

typedef struct edgeNode {
	int iVex, jVex;
	struct edgeNode* iLink, * jLink;
	infoType* info;
}edgeNode;

typedef struct vertexNode {
	int number;
	vertexType data;
	int vi;
	edgeNode* firstEdge;
}vertexNode;

typedef struct {
	vertexNode adjmultiList[vexnumber];
	int vexNum, edgeNum;
}AMLGraph;

//逐边插入；
void insertEdgeAction(AMLGraph& G, int index1, int index2) {
	edgeNode* p, * q;
	edgeNode* edge = new edgeNode[1];

	edge->iVex = index1;
	edge->jVex = index2;

	p = G.adjmultiList[index1].firstEdge;
	if (!p)
	{
		G.adjmultiList[index1].firstEdge = edge;
		edge->iLink = NULL;
	}
	else
	{
		G.adjmultiList[index1].firstEdge = edge;
		edge->iLink = p;
	}

	q = G.adjmultiList[index2].firstEdge;
	if (!q) 
	{
		G.adjmultiList[index2].firstEdge = edge;
		edge->jLink = NULL;
	}
	else
	{
		G.adjmultiList[index2].firstEdge = edge;
		edge->jLink = q;
	}
}
//建图；
void createAMLGraph(AMLGraph& G, int vexNum) {
	G.vexNum = vexNum;
	G.edgeNum = 0;
	for (int i = 0; i < G.vexNum; i++) {
		G.adjmultiList[i].number = i;
		G.adjmultiList[i].data = symbol[i];
		G.adjmultiList[i].vi = 0;
		G.adjmultiList[i].firstEdge = NULL;
	}
	for (int j = 0; j < size(first); j++)
	{
		if (first[j] < 0 || second[j] < 0)
			exit(-1);
		insertEdgeAction(G, first[j], second[j]);
	}
}

void printAMLGraph(AMLGraph& G) {
	for (int i = 0; i < G.vexNum; i++) {
		cout << i << " " << G.adjmultiList[i].data;
		edgeNode* edge = G.adjmultiList[i].firstEdge;

		while (edge) {
			cout << "-->|" << edge->iVex << "|" << edge->jVex << "|";
			if (edge->iVex == i) {
				edge = edge->iLink;
			}
			else {
				edge = edge->jLink;
			}
		}
		cout << "-->NULL" << endl;
	}
}

void BFS(AMLGraph& G, int start)
{
	queue<vertexNode> q;
	q.push(G.adjmultiList[start]);
	while (q.size() > 0)
	{
		vertexNode u = q.front();
		q.pop();
		//访问节点；
		cout << u.data;
		G.adjmultiList[u.number].vi = 1;

		edgeNode* p = u.firstEdge;
		while (p != NULL)
			if (u.data == G.adjmultiList[p->iVex].data && G.adjmultiList[p->jVex].vi == 0)
			{
				q.push(G.adjmultiList[p->jVex]);
				p = p->iLink;
			}
			else 
				p = NULL;
	}
}

int main(int argc)
{
	AMLGraph G;
	createAMLGraph(G, vexnumber);
	cout << "邻接多重表表示：" << endl;
	printAMLGraph(G);

	cout << endl << "BFS广度优先遍历：";
	//BFS(G, 0);
	for (int i = 0; i < vexnumber; i++)
		if (G.adjmultiList[i].vi == 0)
			BFS(G, i);
	cout << endl;

}

