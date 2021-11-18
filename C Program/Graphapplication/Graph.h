#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define OK 1
#define False 0
#define Max 10000

typedef struct {
	char vexs[100];
	int arcs[100][100];
	int vexnum,arcnum;
}Graph;
typedef struct QNode{
	int element;
	struct QNode *next;
}Qnode;
typedef struct Queue{
	QNode *front,*rear;
}Queue;

int LocateVex(Graph G,char v){
	int i;
	for(i=0;i<G.vexnum;++i)
		if(G.vexs[i]==v)
			return i;
	return -1;
}
int CreateUDN(Graph &G){
	int i,j,k;
	char v1,v2;
	printf("Please enter the number of vertex:");
	scanf("%d",&G.vexnum);
	printf("Please enter the number of arc:");
	scanf("%d",&G.arcnum);
	printf("Please enter the detail of vertex in order:");
	getchar();
	for(i=0;i<G.vexnum;i++) scanf("%c",&G.vexs[i]);
	for(i=0;i<G.vexnum;i++)
		for(j=0;j<G.vexnum;j++) G.arcs[i][j]=0;
	//printf("Please enter the edge vertex pair:");
	//getchar();
	//scanf("%c%c",&v1,&v2);
	for(k=0;k<G.arcnum;++k){
		printf("Please enter the edge vertex pair:");
	    getchar();
	    scanf("%c%c",&v1,&v2);
		i=LocateVex(G,v1);
		j=LocateVex(G,v2);
		G.arcs[i][j]=1;
		G.arcs[j][i]=1;
		//printf("Please enter the edge vertex pair:");
	    //scanf("%c %c",&v1,&v2);
	}
	return OK;
}
int FirstAdjVex(Graph G,int v){
	int i;
	if(v<0 || v>=G.vexnum) return -1;
	for(i=0;i<G.vexnum;i++){
		if(G.arcs[v][i]) return i;
	}
	return -1;
}
int NextAdjVex(Graph G,int v,int w){
	if(v<0 || v>=G.vexnum || w>=G.vexnum) return -1;
	for(w=w+1;w<G.vexnum;w++)
		if(G.arcs[v][w]) return w;
	return -1;
}

int visited[100];//访问数组

void DFS(Graph G,int v){
	int w;
	printf("%c",G.vexs[v]);
	visited[v]=OK;
	for(w=FirstAdjVex(G,v);w!=-1;w=NextAdjVex(G,v,w))
		if(!visited[w]) DFS(G,w);
}

void DFSTraverse(Graph G){
	int i,v;
	for(i=0;i<G.vexnum;i++) visited[i]=0;
	for(v=0;v<G.vexnum;v++)
		if(!visited[v]) DFS(G,v);
}
int InsertArc(Graph &G){
	char v1,v2;
	int re=1,i,j;
	printf("Please enter the edge vertex pair that you want insert the arc:");
	getchar();
	scanf("%c%c",&v1,&v2);
	while(re){
		//printf("Please enter the edge vertex pair that you want insert the arc:");
		//scanf("%c%c",&v1,&v2);
		i=LocateVex(G,v1);
		j=LocateVex(G,v2);
		G.arcs[i][j]=1;
		G.arcs[j][i]=1;
		printf("Do you want to enter another Arc?(Yes is 1;No is 0):");
		getchar();
		scanf("%d",&re);
		if(!re) break;
		printf("Please enter the edge vertex pair that you want insert the arc:");
		getchar();
		scanf("%c%c",&v1,&v2);
	}
	return OK;
}
int DeleteArc(Graph &G){
	char v1,v2;
	int re=1,i,j;
	printf("Please enter the edge vertex pair that you want delete the arc:");
	getchar();
	scanf("%c%c",&v1,&v2);
	while(re){
		//printf("Please enter the edge vertex pair that you want delete the arc:");
		//scanf("%c%c",&v1,&v2);
		i=LocateVex(G,v1);
		j=LocateVex(G,v2);
		G.arcs[i][j]=0;
		G.arcs[j][i]=0;
		printf("Do you want to delete another Arc?(Yes is 1;No is 0):");
		getchar();
		scanf("%d",&re);
		if(!re) break;
		printf("Please enter the edge vertex pair that you want delete the arc:");
		getchar();
		scanf("%c%c",&v1,&v2);
	}
	return OK;
}
int DegreeVex(Graph G,char v){
	int i,j,degree=0;
	i=LocateVex(G,v);
	for(j=0;j<G.vexnum;j++)
		if(G.arcs[i][j]) degree++;
	printf("The degree of the vertex '%c' is %d",v,degree);
	return OK;
}
//部分所需队列函数.
int InitQueue(Queue &Q){
	Q.front=Q.rear=(QNode *)malloc(sizeof(QNode));
	if(!Q.front) return False;
	Q.front->next=NULL;
	return OK;
}
int QueueEmpty(Queue Q){
	if(Q.front==Q.rear) return OK;
	else return False;
}
int EnQueue(Queue &Q,int e){
	QNode *p=NULL;
	p=(QNode *)malloc(sizeof(QNode));
	if(!p) return False;
	p->element=e;
	Q.rear->next=p;
	Q.rear=p;
	return OK;
}
int DeQueue(Queue &Q,int &e){
	QNode *p;
	if(Q.front==Q.rear) return False;
	p=Q.front->next;
	e=p->element;
	Q.front->next=p->next;
	if(Q.rear==p) Q.rear=Q.front;
	free(p);
	return OK;
}
//OVER.
void BFSTraverse(Graph G){
	Queue Q;
	int v,w,u;
	for(v=0;v<G.vexnum;v++) visited[v]=0;
	InitQueue(Q);
	for(v=0;v<G.vexnum;v++){
		if(!visited[v]){
			printf("%c",G.vexs[v]);
			visited[v]=OK;
			EnQueue(Q,v);
			while(!QueueEmpty(Q)){
				DeQueue(Q,u);
				w=FirstAdjVex(G,u);
				while(w!=-1){
					if(!visited[w]){
						printf("%c",G.vexs[w]);
						visited[w]=OK;
						EnQueue(Q,w);
					}
					w=NextAdjVex(G,u,w);
				}
			}
		}
	}
}
int MinPath(Graph G,char v,char x){
	int i,m,n,u,w;
	Queue Q;
	m=LocateVex(G,v);
	n=LocateVex(G,x);
	for(i=0;i<G.vexnum;i++) visited[i]=-1;
	InitQueue(Q);
	visited[m]=0;
	EnQueue(Q,m);
	while(!QueueEmpty(Q)){
		DeQueue(Q,u);
		for(w=FirstAdjVex(G,u);w!=-1;w=NextAdjVex(G,u,w)){
			if(visited[w]!=-1) continue;
			visited[w]=visited[u]+1;
			if(w==n){
				printf("The min path is %d",visited[n]);
				return OK;
			}
			EnQueue(Q,w);
		}
	}
}
int CreateUDN_1(Graph &G){
	int i,j,k,weight;
	char v1,v2;
	printf("Please enter the number of vertex:");
	scanf("%d",&G.vexnum);
	printf("Please enter the number of arc:");
	scanf("%d",&G.arcnum);
	printf("Please enter the detail of vertex in order:");
	getchar();
	for(i=0;i<G.vexnum;i++) scanf("%c",&G.vexs[i]);
	for(i=0;i<G.vexnum;i++)
		for(j=0;j<G.vexnum;j++) G.arcs[i][j]=Max;
	//printf("Please enter the edge vertex pair:");
	//getchar();
	//scanf("%c%c",&v1,&v2);
	for(k=0;k<G.arcnum;++k){
		printf("Please enter the edge vertex pair and the weight of the Arc:");
	    getchar();
	    scanf("%c%c %d",&v1,&v2,&weight);
		i=LocateVex(G,v1);
		j=LocateVex(G,v2);
		G.arcs[i][j]=weight;
		G.arcs[j][i]=weight;
		//printf("Please enter the edge vertex pair:");
	    //scanf("%c %c",&v1,&v2);
	}
	return OK;
}
void Dijkstra(Graph G,char v0){
	int S[100]={0},D[100],P[100][100],i,j,k,w,m,min;
	char v;
	m=LocateVex(G,v0);
	for(i=0;i<G.vexnum;i++){
		D[i]=G.arcs[m][i];
		if(D[i]!=Max){
			P[i][0]=m;
			P[i][1]=i;
			P[i][2]=-1;
		}
	}
	S[m]=1;
	D[m]=0;
	for(i=0;i<G.vexnum;i++){
		min=Max;
		for(j=0;j<G.vexnum;j++)
			if(!S[j] && D[j]<min){
				min=D[j];
				k=j;
			}
		S[k]=1;
		for(j=0;j<G.vexnum;j++)
			if(!S[j] && D[k]+G.arcs[k][j]<D[j]){
				D[j]=D[k]+G.arcs[k][j];
				for(w=0;P[k][w]!=-1;w++) P[j][w]=P[k][w];
				P[j][w]=j;
				P[j][w+1]=-1;
			}
	}
	//for(i=0;i<G.vexnum;i++){
	//	for(j=0;j<=G.vexnum;j++) 
	//		if(P[i][j]>=0) printf("%d",P[i][j]);
	//	printf("\n");
	//}
	//printf("Please enter the out vertex:");
	//getchar();
	//scanf("%c",&v);
	//j=LocateVex(G,v);
	for(j=0;j<G.vexnum;j++){
		for(k=0;k<=G.vexnum;k++) 
			if(P[j][k]>=0) printf("%c",G.vexs[P[j][k]]);
		printf("\n");
	}
	//for(i=0;i<G.vexnum;i++) printf("%d\n",S[i]);
	//for(i=0;i<G.vexnum;i++) printf("%d\n",D[i]);
}
int MinEdge(int a[],int n){
	int i,min=Max,j=i;
	for(i=0;i<n;i++)
		if(a[i]!=0 && a[i]<min){
			min=a[i];
			j=i;
		}
	return j;
}
void Prim(Graph G,char v){
	int lowcost[100]={0},v0,i,j,k;
	char adjvex[100]={'\0'};
	v0=LocateVex(G,v);
	lowcost[v0]=0;
	for(j=0;j<G.vexnum;j++)
		if(j!=v0){
			lowcost[j]=G.arcs[v0][j];
			adjvex[j]=v;
		}
	for(i=1;i<G.vexnum;i++){
		k=MinEdge(lowcost,G.vexnum);
		printf("(%c,%c)\n",adjvex[k],G.vexs[k]);
		lowcost[k]=0;
		for(j=0;j<G.vexnum;j++)
			if(G.arcs[k][j]<lowcost[j]){
				adjvex[j]=G.vexs[k];
				lowcost[j]=G.arcs[k][j];
			}
	}
}



int ShortestPath(Graph G){//每队顶点间的最小路径
	int D[100][100],P[100][100][100],v,u,w,i,j;
	char v1,v2;
	for(v=0;v<G.vexnum;++v)
		for(w=0;w<G.vexnum;++w){
			D[v][w]=G.arcs[v][w];
			for(u=0;u<G.vexnum;++u) P[v][w][u]=-1;
			if(D[v][w]<Max){
				P[v][w][v]=OK;
				P[v][w][w]=OK;
			}
		}
	for(u=0;u<G.vexnum;++u)
		for(v=0;v<G.vexnum;++v)
			for(w=0;w<G.vexnum;++w)
				if(D[v][u]+D[u][w] < D[v][w]){
					D[v][w]=D[v][u]+D[u][w];
					for(i=0;i<G.vexnum;++i)
						P[v][w][i]=P[v][u][i] || P[u][w][i];
				}
	printf("Please enter the two vertex:");
	getchar();
	scanf("%c%c",&v1,&v2);
	i=LocateVex(G,v1);
	j=LocateVex(G,v2);
	for(v=0;v<G.vexnum;v++)
		printf("%c\n",G.vexs[P[i][j][v]]);
	return OK;
}