#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define Max 100
#define OK 1
#define False 0

typedef struct BiTree{
	char data;
	struct BiTree *lchild,*rchild;
}BiTree;
typedef struct QNode{
	char element;
	struct QNode *next;
}Qnode;
typedef struct Queue{
	QNode *front,*rear;
}Queue;

BiTree *InilBiTree(BiTree *T,char **a){
	//printf("%c\n", **a);
	if(**a == '#' || **a == '\0'){
		return NULL;
	}
	else{
		T=(BiTree *)malloc(sizeof(BiTree));
		T->lchild=NULL;
		T->rchild=NULL;
		T->data=**a;
	}
	++*a;
	//printf("%c",*a);
	T->lchild=InilBiTree(T->lchild,a);
	//printf("a");
	++*a;
	//printf("%c",*a);
	T->rchild=InilBiTree(T->rchild,a);
	return T;
}
int DestroyBiTree(BiTree *T){
	if(!T) return OK;
	DestroyBiTree(T->lchild);
	DestroyBiTree(T->rchild);
	free(T);
	return OK;
}
int PreOrder(BiTree *T){
	if(!T) return OK;
	printf("%c",T->data);
	PreOrder(T->lchild);
	PreOrder(T->rchild);
	return OK;
}
int InOrder(BiTree *T){
	if(!T) return OK;
	InOrder(T->lchild);
	printf("%c",T->data);
	InOrder(T->rchild);
	return OK;
}
int PostOrder(BiTree *T){
	if(!T) return OK;
	PostOrder(T->lchild);
	PostOrder(T->rchild);
	printf("%c",T->data);
	return OK;
}
void Nodes(BiTree *T,int &c){
	if(!T) return;
	c++;
	Nodes(T->lchild,c);
	Nodes(T->rchild,c);
}
void Leaves(BiTree *T,int &c){
	if(!T) return;
	if(T->lchild==NULL && T->rchild==NULL) c++;
	Leaves(T->lchild,c);
	Leaves(T->rchild,c);
}
int Depth(BiTree *T){
	int l,r;
	if(!T) return 0;
	l=Depth(T->lchild);
	r=Depth(T->rchild);
	return (l>r? l+1 : r+1);
}
int Search(BiTree *T,char x,char &parent,char &lch,char &rch,char &brother){
	if(!T) return OK;
	if(T->lchild && T->lchild->data==x){
		parent=T->data;
		if(T->rchild) brother=T->rchild->data;
		if(T->lchild->lchild) lch=T->lchild->lchild->data;
		if(T->lchild->rchild) rch=T->lchild->rchild->data;
		return OK;
	}
	if(T->rchild && T->rchild->data==x){
		parent=T->data;
		if(T->lchild) brother=T->lchild->data;
		if(T->rchild->lchild) lch=T->rchild->lchild->data;
		if(T->rchild->rchild) rch=T->rchild->rchild->data;
		return OK;
	} else{
		Search(T->lchild,x,parent,lch,rch,brother);
		Search(T->rchild,x,parent,lch,rch,brother);
		return OK;
	}
}
BiTree *CopyTree(BiTree *T){
	BiTree *newlp=NULL,*newrp=NULL,*newnode=NULL;
	if(!T) return NULL;
	newlp=(BiTree *)malloc(sizeof(BiTree));
	newrp=(BiTree *)malloc(sizeof(BiTree));
	newnode=(BiTree *)malloc(sizeof(BiTree));
	newlp=CopyTree(T->lchild);
	newrp=CopyTree(T->rchild);
	newnode->data=T->data;
	newnode->lchild=newlp;
	newnode->rchild=newrp;
	return newnode;
}
BiTree *Location(BiTree *T,char x){//定位
	BiTree *k=NULL;
	if(!T) return NULL;
	if(T->data==x) return T;
	k=Location(T->lchild,x);
	if(!k) k=Location(T->rchild,x);
	else return k;
}
int Delete(BiTree *T,char x,int direction){
	BiTree *p=NULL;
	p=Location(T,x);
	if(!p){
		printf("The node that needs to be deleted doesn't exist.\n");
		return OK;
	}
	switch(direction){
	case 1:
		if(!p->lchild){
			printf("The node that needs to be deleted doesn't exist.\n");
			return OK;
		}
		DestroyBiTree(p->lchild);
		p->lchild=NULL;
		printf("Deleted successfully.\n");
		printf("The preordering sequence of the deleted BiTree is");
		PreOrder(T);
		return OK;
	case 2:
		if(!p->rchild){
			printf("The node that needs to be deleted doesn't exist.\n");
			return OK;
		}
		DestroyBiTree(p->rchild);
		p->rchild=NULL;
		printf("Deleted successfully.\n");
		printf("The preordering sequence of the deleted BiTree is");
		PreOrder(T);
		return OK;
	}
}
int Insert(BiTree *T,BiTree *insert,char x,int direction){
	BiTree *p=NULL;
	p=Location(T,x);
	if(!insert){
		printf("The subtree to be inserted is NULL.\n");
		return OK;
	}
	if(!p){
		printf("The node that needs to insert doesn't exist.\n");
		return OK;
	}
	switch(direction){
	case 1:
		if(p->lchild){
			printf("The position to be inserted isn't empty.\n");
			return OK;
		}
		p->lchild=insert;
		printf("Inserted successfully.\n");
		printf("The preordering sequence of the inserted BiTree is");
		PreOrder(T);
		return OK;
	case 2:
		if(p->rchild){
			printf("The position to be inserted isn't empty.\n");
			return OK;
		}
		p->rchild=insert;
		printf("Inserted successfully.\n");
		printf("The preordering sequence of the inserted BiTree is");
		PreOrder(T);
		return OK;
	}
}
int Exchange(BiTree *T,char x){
	BiTree *p=NULL,*q=NULL;
	p=Location(T,x);
	if(!p){
		printf("The node that needs to exchange doesn't exist.\n");
		return OK;
	}
	if(p->lchild==NULL && p->rchild==NULL) return OK;
	if(p->lchild==NULL && p->rchild){
		p->lchild=p->rchild;
		DestroyBiTree(p->rchild);
		printf("Exchanged successfully.\n");
		printf("The preordering sequence of the exchanged BiTree is");
		PreOrder(T);
		return OK;
	}
	if(p->lchild && p->rchild==NULL){
		p->rchild=p->lchild;
		DestroyBiTree(p->lchild);
		printf("Exchanged successfully.\n");
		printf("The preordering sequence of the exchanged BiTree is");
		PreOrder(T);
		return OK;
	}
	q=p->lchild;
	p->lchild=p->rchild;
	p->rchild=q;
	printf("Exchanged successfully.\n");
	printf("The preordering sequence of the exchanged BiTree is");
	PreOrder(T);
	return OK;
}
int Judge_BiSortTree(BiTree *T){
	if(!T || (!T->lchild && !T->rchild)) return OK;
	if(T->lchild->data >= T->data || T->rchild->data < T->data) return False;
	if(!Judge_BiSortTree(T->lchild)) return False;
	if(!Judge_BiSortTree(T->rchild)) return False;
	return OK;
}
int Judge_AVL_BiTree(BiTree *T){
	int result;
	if(!T || (!T->lchild && !T->rchild)) return OK;
	result=Depth(T->lchild) - Depth(T->rchild);
	switch(result){
	case -1:
	case 0:
	case 1:
		if(Judge_AVL_BiTree(T->lchild) && Judge_AVL_BiTree(T->rchild)) return OK;
		else return False;
	default :return False;
	}
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
int EnQueue(Queue &Q,char e){
	QNode *p=NULL;
	p=(QNode *)malloc(sizeof(QNode));
	if(!p) return False;
	p->element=e;
	Q.rear->next=p;
	Q.rear=p;
	return OK;
}
int DeQueue(Queue &Q,char &e){
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
void LayerTaversal(BiTree *T){
	char p;
	BiTree *local=NULL;
	Queue Q;
	InitQueue(Q);
	if(T) EnQueue(Q,T->data);
	while(!QueueEmpty(Q)){
		DeQueue(Q,p);
		local=Location(T,p);
		printf("%c",local->data);
		if(local->lchild) EnQueue(Q,local->lchild->data);
		if(local->rchild) EnQueue(Q,local->rchild->data);
	}
}
int Judge_Complete_BiTree(BiTree *T){
	char p,o='#';
	BiTree *local=T;
	Queue Q;
	InitQueue(Q);
	if(!T) return OK;
	EnQueue(Q,T->data);
	while(local){
		DeQueue(Q,p);
		if(p=='#') break;
		local=Location(T,p);
		//if(!local->lchild || !local->rchild) return False;
		if(local){
			if(local->lchild) EnQueue(Q,local->lchild->data);
			else EnQueue(Q,o);
			if(local->rchild) EnQueue(Q,local->rchild->data);
			else EnQueue(Q,o);
		}
	}
	
	while(!QueueEmpty(Q) && p=='#'){
		DeQueue(Q,p);
		//if(p=='#') return False;
	}
	if(p!='#') return False;
	else return OK;
	
}
