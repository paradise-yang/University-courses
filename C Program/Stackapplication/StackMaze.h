#include<stdio.h>
#include<stdlib.h>

#define OK 1
#define False 0
typedef struct Pos {//迷宫单个结点定义
	int row,column;
}Pos;
typedef struct SElemType {//迷宫路径
	int ord;
	Pos seat;
	int direction;
}SElemType;
typedef char MazeType[10][10];

typedef struct MStack {
	SElemType *base;
	SElemType *top;
	int stacksize;
}MStack;

int InitStack(MStack &S){
	S.base = (SElemType *)malloc(Stack_init_size*sizeof(SElemType));
	if(!S.base) return False;
	S.top=S.base;
	S.stacksize=Stack_init_size;
	return OK;
}
int StackEmpty(MStack S){
	if(S.top==S.base) return OK;
	else return False;
}
//int GetTop(MStack S,SElemType *e){
//	if(S.top==S.base){
//		printf("栈空\n");
//		return False;
//	}
//	else{
//		e=*(S.top-1);
//		return OK;
//	}
//}
int Push(MStack &S,SElemType e){
	if(S.top-S.base>=S.stacksize){
		S.base=(SElemType *)realloc(S.base,(S.stacksize+Stackincrement)*sizeof(SElemType));
		if(!S.base) exit(0);
		S.top=S.base+S.stacksize;
		S.stacksize+=Stackincrement;
	}
	*S.top++ = e;
	return OK;
}
int Pop(MStack &S,SElemType &e){
	if(S.top==S.base) return False;
	e=*--S.top;
	return OK;
}
int DestroyStack(MStack &S){
	if(StackEmpty(S)) return False;
	free(S.top);
	free(S.base);
	return OK;
}

//迷宫函数

void FootPrint(MazeType &maze,Pos &p){//留下足迹
	maze[p.row][p.column]='*';
}
void MarkPrint(MazeType &maze,Pos &p){//不能通行
	maze[p.row][p.column]='X';
}
int Pass(MazeType &maze,Pos &p){//标记通过
	return (maze[p.row][p.column]=='0');
}
Pos NextPos(Pos p,int d){//下
	Pos n;
	n.column=p.column+((d==1)?1:(d==3)?-1:0);
	n.row=p.row+((d==2)?1:(d==4)?-1:0);
	return n;
}
int MazePath(MazeType &maze,Pos start,Pos end){//迷宫求解
	Pos curpos;
	MStack S;
	SElemType e;
	int curstep=1;
	curpos.row=start.row;
	curpos.column=start.column;
	InitStack(S);
	do{
		if(Pass(maze,curpos)==1){
			FootPrint(maze,curpos);
			e.ord=curstep;
			e.seat.row=curpos.row;
			e.seat.column=curpos.column;
			e.direction=1;
			Push(S,e);
			if(curpos.row==end.row && curpos.column==end.column) return OK;
			curpos=NextPos(curpos,1);
			curstep++;
		}
		else {
			if(!StackEmpty(S)){
				Pop(S,e);
				while(e.direction==4 && !StackEmpty(S)){
					MarkPrint(maze,e.seat);
					Pop(S,e);
				}
				if(e.direction<4){
					e.direction++;
					Push(S,e);
					curpos=NextPos(e.seat,e.direction);
				}
			}
		}
	}while(!StackEmpty(S));
	printf("There is no path.\n");
	return False;
}
