#include<stdio.h>
#include<stdlib.h>
#include<math.h>
typedef struct text{
	float coeff;          //系数
	int index;            //指数
    struct text *next;
}Polyn;
Polyn *InsertPolyn(Polyn *polyn,Polyn *e){                  //插入函数
	Polyn *p,*q;
	p=polyn;
	q=polyn;
	if(polyn==NULL){
		polyn=e;
		e->next=NULL;
		return  polyn;
	}
	if(e->index > p->index){
		e->next=polyn;
		polyn=e;
		return polyn;
	}
	while(p!=NULL && e->index < p->index){
		q=p;
		p=p->next;
	}
	if(p==NULL){
		q->next=e;
		e->next=NULL;
		return polyn;
	}
	if(p->index==e->index){
		p->coeff=p->coeff + e->coeff;
		if(p->coeff==0){
			q->next=p->next;
			free(p);
			return polyn;
		}
		return polyn;
	}
	q->next=e;
	e->next=p;
	return polyn;
}
Polyn *CreatPolyn(void){                                    //创建多项式
	Polyn *polyn=NULL,*p;
	float xs=0;
	int zs=0;
	printf("\n请输入首项系数：");
	scanf("%f",&xs);
	if(xs==0) return NULL;
	while(xs!=0){
		p=(Polyn *)malloc(sizeof(Polyn));
		p->coeff=xs;
		printf("请输入该项指数:");
		scanf("%d",&p->index);
		polyn=InsertPolyn(polyn,p);
		printf("请输入下一项系数（输入0结束）:");
		scanf("%f",&xs);
	}
	printf("多项式建立完成。\n");
	return polyn;
}
int PrintPolyn(Polyn *polyn){                               //打印多项式
	Polyn *p;
	p=polyn;
	if(!p){
		printf("多项式为空！");
		return(0);
	}
	while(p!=NULL){
		if(p->coeff >0) printf("+%.2f",p->coeff);
		else printf("%.2f",p->coeff);
		if(p->index !=0) printf("x^%d",p->index);
		else printf("\n");
		p=p->next;
	}
	printf("\n");
	return(0);
}
void DestroyPolyn(Polyn *polyn){                            //销毁多项式
	Polyn *p,*q;
	q=polyn;
	p=q->next;
	while(p!=NULL){
		free(q);
		q=p;
		p=p->next;
	}
	free(q);
}
void AddPolyn(Polyn *polyn1,Polyn *polyn2){                 //多项式加法
	Polyn *p,*q,*tempolyn=NULL,*rear=NULL,*t;
	p=polyn1;
	q=polyn2;
	while(p!=NULL && q!=NULL){
		t=(Polyn *)malloc(sizeof(Polyn));
		if(p->index > q->index){
			t->coeff=p->coeff;
			t->index=p->index;
			p=p->next;
		}
		else if(p->index < q->index){
			t->coeff=q->coeff;
			t->index=q->index;
			q=q->next;
		}
		else {
			t->coeff=p->coeff + q->coeff;
			if(t->coeff!=0){
				t->index=p->index;
				p=p->next;
				q=q->next;
			}
			else{
				free(t);
				p=p->next;
				q=q->next;
				continue;
			}
		}
		if(tempolyn==NULL) tempolyn=t;
		else rear->next=t;
		rear=t;
	}
	if(p!=NULL){
		while(p!=NULL){
			t=(Polyn *)malloc(sizeof(Polyn));
			t->coeff=p->coeff;
			t->index=p->index;
			if(tempolyn==NULL) tempolyn=t;
			else rear->next=t;
			rear=t;
			p=p->next;
		}
		if(rear!=NULL) rear->next=NULL;
	}
	if(q!=NULL){
		while(q!=NULL){
			t=(Polyn *)malloc(sizeof(Polyn));
			t->coeff=q->coeff;
			t->index=q->index;
			if(tempolyn==NULL) tempolyn=t;
			else rear->next=t;
			rear=t;
			q=q->next;
		}
		if(rear!=NULL) rear->next=NULL;
	}
	if(rear!=NULL) rear->next=NULL;
	printf("结果如下：\n");
	PrintPolyn(tempolyn);
	DestroyPolyn(tempolyn);
}
Polyn *SubtractPolyn(Polyn *polyn1,Polyn *polyn2){          //多项式减法
	Polyn *p,*q,*tempolyn=NULL,*rear=NULL,*t;
	p=polyn1;
	q=polyn2;
	while(p!=NULL && q!=NULL){
		t=(Polyn *)malloc(sizeof(Polyn));
		if(p->index > q->index){
			t->coeff=p->coeff;
			t->index=p->index;
			p=p->next;
		}
		else if(p->index < q->index){
			t->coeff=-(q->coeff);
			t->index=q->index;
			q=q->next;
		}
		else {
			t->coeff=p->coeff - q->coeff;
			if(t->coeff!=0){
				t->index=p->index;
				p=p->next;
				q=q->next;
			}
			else{
				free(t);
				p=p->next;
				q=q->next;
				continue;
			}
		}
		if(tempolyn==NULL) tempolyn=t;
		else rear->next=t;
		rear=t;
	}
	if(p!=NULL){
		while(p!=NULL){
			t=(Polyn *)malloc(sizeof(Polyn));
			t->coeff=p->coeff;
			t->index=p->index;
			if(tempolyn==NULL) tempolyn=t;
			else rear->next=t;
			rear=t;
			p=p->next;
		}
		if(rear!=NULL) rear->next=NULL;
	}
	if(q!=NULL){
		while(q!=NULL){
			t=(Polyn *)malloc(sizeof(Polyn));
			t->coeff=-(q->coeff);
			t->index=q->index;
			if(tempolyn==NULL) tempolyn=t;
			else rear->next=t;
			rear=t;
			q=q->next;
		}
		if(rear!=NULL) rear->next=NULL;
	}
	if(tempolyn==NULL) return NULL;
	if(rear!=NULL) rear->next=NULL;
	return tempolyn;
}
double EvaluePolyn(Polyn *polyn,float x){                   //多项式求值
	Polyn *p;
	double sum=0;
	p=polyn;
	while(p!=NULL){
		sum=sum+(p->coeff)*pow(x,p->index);
		p=p->next;
	}
	return sum;
}
void DifferentiatePolyn(Polyn *polyn,int n){                //求微分
	Polyn *p,*t,*tempolyn=NULL,*rear;
	int i;
	p=polyn;
	while(p!=NULL){
		t=(Polyn *)malloc(sizeof(Polyn));
		
		if(p->index >= n){
			t->coeff=p->coeff;
			for(i=0;i<n;i++) t->coeff=(t->coeff)*((p->index)-i);
			t->index=p->index-n;
			p=p->next;
		}
		else{
			p=p->next;
			free(t);
			continue;
		}
		
		if(tempolyn==NULL) tempolyn=t;
		else rear->next=t;
		rear=t;
	}
	if(rear!=NULL) rear->next=NULL;
	PrintPolyn(tempolyn);
	DestroyPolyn(tempolyn);
}
double Definite_integralPolyn(Polyn *polyn,float a,float b){//求定积分
	Polyn *p,*t,*tempolyn=NULL,*rear;
	double uplimit=0,lowlimit=0,result=0;
	p=polyn;
	while(p!=NULL){
		t=(Polyn *)malloc(sizeof(Polyn));
		t->index=p->index+1;
		t->coeff=p->coeff/(t->index);
		p=p->next;
		if(tempolyn==NULL) tempolyn=t;
		else rear->next=t;
		rear=t;
	}
	if(rear!=NULL) rear->next=NULL;
	uplimit=EvaluePolyn(tempolyn,a);
	lowlimit=EvaluePolyn(tempolyn,b);
	result=uplimit - lowlimit;
	return result;
}
Polyn *MultiplyPolyn(Polyn *polyn1,Polyn *polyn2){          //多项式乘法
	Polyn *p,*q,*k,*tempolyn=NULL;
	p=polyn1;
	q=polyn2;
	while(p!=NULL){
		while(q!=NULL){
			k=(Polyn *)malloc(sizeof(Polyn));
			if(k==NULL){
				printf("ERROR");
				return (0);
			}
			k->index=p->index + q->index;
			k->coeff=p->coeff * q->coeff;
			tempolyn=InsertPolyn(tempolyn,k);
			q=q->next;
		}
		q=polyn2;
		p=p->next;
	}
	return tempolyn;
}
Polyn *CopyPolyn(Polyn *polyn){                             //复制函数
	Polyn *tempolyn=NULL,*rear,*p,*t;
	p=polyn;
	while(p!=NULL){
		t=(Polyn *)malloc(sizeof(Polyn));
		t->coeff=p->coeff;
		t->index=p->index;
		if(tempolyn==NULL) tempolyn=t;
		else rear->next=t;
		rear=t;
		p=p->next;
	}
	if(rear!=NULL) rear->next=NULL;
	return tempolyn;
}
int DivisionPolyn(Polyn *polyn1,Polyn *polyn2){             //多项式除法
	Polyn *p,*t,*rear;
	Polyn *quotpolyn=NULL,*remainpolyn=NULL,*tempolyn=NULL;
	remainpolyn=CopyPolyn(polyn1);
	p=polyn1;
	if(polyn1->index < polyn2->index){
		printf("商多项式：0");
		printf("余多项式：");
		PrintPolyn(remainpolyn);
		return(0);
	}
	while(p->index >= polyn2->index){
		t=(Polyn *)malloc(sizeof(Polyn));
		if(t==NULL){
				printf("ERROR");
				return (0);
			}
		t->coeff=p->coeff / polyn2->coeff;
		t->index=p->index - polyn2->index;
		t->next=NULL;
		if(quotpolyn==NULL) quotpolyn=t;
		else rear->next=t;
		rear=t;
        remainpolyn=SubtractPolyn(remainpolyn,MultiplyPolyn(t,polyn2));
		p=remainpolyn;
		if(p==NULL) break;
	}
	if(!rear) rear->next=NULL;
	printf("商多项式：");
	PrintPolyn(quotpolyn);
	printf("余多项式：");
	if(p==NULL) printf("0.00\n");
	else PrintPolyn(remainpolyn);
}
