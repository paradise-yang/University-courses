#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define Stack_init_size 100
#define Stackincrement  10
#define OK 1
#define False 0

struct Stack {
	char *base;
	char *top;
	int stacksize;
};

int InitStack(Stack &S){
	S.base = (char *)malloc(Stack_init_size*sizeof(char));
	if(!S.base) return False;
	S.top=S.base;
	S.stacksize=Stack_init_size;
	return OK;
}
int StackEmpty(Stack S){
	if(S.top==S.base) return OK;
	else return False;
}

char GetTop(Stack S,char &e){
	if(S.top==S.base) return False;
	e=*(S.top-1);
	return e;
}
int Push(Stack &S,char e){
	if(S.top-S.base>=S.stacksize){
		S.base=(char *)realloc(S.base,(S.stacksize+Stackincrement)*sizeof(char));
		if(!S.base) exit(0);
		S.top=S.base+S.stacksize;
		S.stacksize+=Stackincrement;
	}
	*S.top++ = e;
	return OK;
}
int Pop(Stack &S,char &e){
	if(S.top==S.base) return False;
	e=*--S.top;
	return OK;
}
int TraverseStack(Stack S){
	if(StackEmpty(S)){
		printf("’ªø’\n");
		return False;
	}
	while(S.base==S.top) printf("%c",*S.base++);
	return OK;
}
int DestroyStack(Stack &S){
	if(StackEmpty(S)) return False;
	free(S.top);
	free(S.base);
	return OK;
}
int BracketsMatching(char *c){
	Stack S;
	char e;
	int state=1;
	InitStack(S);//≥ı º’ª
	e=*c++;
	while(e!='\0' && state){
		switch(e){
		case '(':
		case '[':
		case '{':
			Push(S,e);
			break;
		case ')':
			GetTop(S,e);
			if(!StackEmpty(S) && e=='(') Pop(S,e);
			else state=0;
			break;
		case ']':
			GetTop(S,e);
			if(!StackEmpty(S) && e=='[') Pop(S,e);
		    else state=0;
			break;
		case '}':
			GetTop(S,e);
			if(!StackEmpty(S) && e=='{') Pop(S,e);
			else state=0;
			break;
		  }
		e=*c++;

	}
	if(state==1 && StackEmpty(S)) printf("The brackets are correct!\n");
	else printf("The brackets are false!\n");
	return OK;
}
