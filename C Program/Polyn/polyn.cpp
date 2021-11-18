#include"Polyn.h"

void main(){
	int request,count;
	float x,a,b;
	Polyn *polyn1=NULL,*polyn2=NULL,*tempolyn=NULL;
	printf("*****����ʽ��������*****");
	while(1){
		printf("\n0.Exit.\n");
		printf("1.Creat a polynomial.\n");
		printf("2.Print a polynomial.\n");
		printf("3.Polynomial summation.\n");
		printf("4.Polynomial differentiation.\n");
		printf("5.Evaluate value.\n");
		printf("6.Differentiate a polynomial.\n");
		printf("7.Take the definite integral of a polynomial.\n");
		printf("8.Polynomial multiplication.\n");
		printf("9.Polynomial division.\n");
		printf("10.Destroy a polynomial.\n");
		printf("Please input your request:");
		scanf("%d",&request);
		if(request==0) break;
		switch(request){
		case 1:
			polyn1=CreatPolyn();
			polyn2=CreatPolyn();
			break;
		case 2:
			printf("polyn1=");
			PrintPolyn(polyn1);
			printf("polyn2=");
			PrintPolyn(polyn2);
			break;
		case 3:
			AddPolyn(polyn1,polyn2);
			break;
		case 4:
			tempolyn=SubtractPolyn(polyn1,polyn2);
		    printf("������£�");
			PrintPolyn(tempolyn);
			DestroyPolyn(tempolyn);
			break;
		case 5:
		    printf("�������Ա���x=");
		    scanf("%f",&x);
		    printf("polyn1(%5.2f)=%lf\n",x,EvaluePolyn(polyn1,x));
	        printf("polyn2(%5.2f)=%lf\n",x,EvaluePolyn(polyn2,x));
		       break;
		case 6:
	        printf("����������΢�ֽ���n=");
            scanf("%d",&count);
			if(count > polyn1->index || polyn1==NULL) printf("polyn1^(%d)=0.00",count);
			else {
				printf("polyn1^(%d)=",count);
				DifferentiatePolyn(polyn1,count);
			}
			if(count > polyn2->index || polyn2==NULL) printf("polyn2^(%d)=0.00",count);
			else {
				printf("polyn2^(%d)=",count);
		        DifferentiatePolyn(polyn2,count);
			}
		    break;
		case 7:
		    printf("�����붨��������a=");
		    scanf("%f",&a);
		    printf("�����붨��������b=");
		    scanf("%f",&b);
		    printf("polyn1�����ֽ��result=%lf\n",Definite_integralPolyn(polyn1,a,b));
		    printf("polyn2�����ֽ��result=%lf\n",Definite_integralPolyn(polyn2,a,b));
		    break;
		case 8:
		    tempolyn=NULL;
		    tempolyn=MultiplyPolyn(polyn1,polyn2);
		    printf("������£�");
			PrintPolyn(tempolyn);
			DestroyPolyn(tempolyn);
			break;
		case 9:
		    DivisionPolyn(polyn1,polyn2);
		    break;
		case 10:
			DestroyPolyn(polyn1);
			break;
		}
	}
}