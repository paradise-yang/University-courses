#include<stdio.h>
#include<stdlib.h>

#define STACK_INIT_SIZE 100
#define STACKINCREMENT 10

struct SqStack{
	int *base;
	int *top;
	int stacksize;
};

int InitStack(SqStack *S){
	S->base=(int *)malloc(STACK_INIT_SIZE*sizeof(int));
	if(!S->base) return(-1);
	S->top=S->base;
	S->stacksize=STACK_INIT_SIZE;
	return(1);
}
int DestroyStack(SqStack *S){
	free(S->base);
	S->base=NULL;
	S->top=NULL;
	S->stacksize=0;
	if(S->stacksize==0) return(1);
	else return(0);
}
int StackEmpty(SqStack *S){
	if(S->top==S->base) return(1);
	else return(0);
}
int Push(SqStack *S,int e){
	if(S->top-S->base>=S->stacksize){
		S->base=(int *)realloc(S->base,(S->stacksize+STACKINCREMENT)*sizeof(int));
		if(!S->base) return(-1);
		S->top=S->base+S->stacksize;
		S->stacksize=S->stacksize+STACKINCREMENT;
	}
	*S->top=e;
	S->top=S->top+1;
	return(1);
}
int Pop(SqStack *S){
	int e;
	if(S->top==S->base) return(0);
	e=*S->top;
	S->top=S->top-1;
	return(e);
}
int GetTop(SqStack *S){
	int e;
	if(S->top==S->base) return(0);
	else{
		e=*(S->top-1);
		return(e);
	}
}
void TraverseStack(SqStack *S){
	int *c;
	int a[8][9],i,j;
	c=S->base;
	while(S->top>c){
		printf("%4d",*c);
		c++;
	}
	printf("\n");
	c=S->base;
	while(S->top>c){
		for(i=0;i<8;i++){
			for(j=1;j<9;j++){
				if(j==*c) printf("   *");
				else printf("    ");
			}
			c++;
			printf("\n");
		}
	}
	printf("\n");
}
int Judgement1(int a,int b,int c,int d){
	if(b==d||(a+b)==(c+d)||(a-b)==(c-d)) return(0);
	else return(1);
}
int Judgement2(int j,int i,int *a){
	int n,x=1;
	for(n=1;n<i;n++){
		if(Judgement1(i,j,n,a[n])==0){
			x=0;
			break;
		}
	}
	return(x);
}
void Queen(){
	struct SqStack *line;
	int i=2,j,k,x=1,a[9];
	a[0]=a[2]=a[3]=a[4]=a[5]=a[6]=a[7]=a[8]=0;
	line=(struct SqStack *)malloc(sizeof(struct SqStack));
	InitStack(line);
	Push(line,1);
	a[1]=1;
	while(i<=8&&i>0){
		j=a[i]+1;
		while(Judgement2(j,i,a)==0){
			j=j+1;
		}
		if(j>8){
			Pop(line);
			a[i]=0;
			i=i-1;
		}//回到上一行
		else{
			Push(line,j);
			a[i]=j;
			if(i==8){
				TraverseStack(line);
				Pop(line);
			}//满足要求，输出
			else{
				i=i+1;
			}//还差几行
		}
	}
}
