#include"Search.h"

void main(){
	HashTable H;
	int request,element,p,c;
	printf("*****≤È’“*****\n");
	printf("0.Exist.\n");
	printf("1.Creat a HashTable.\n");
	printf("2.Traverse the Hash?Table.\n");
	printf("3.Search.\n");
	printf("4.Calculate ASL.\n");
	printf("5.Delete.\n");
	printf("6.Insert.\n");
	while(1){
		printf("\nPlease enter your request:");
		scanf("%d",&request);
		switch(request){
		case 1:
			InitDsTable(H);
			break;
		case 2:
			Traverse(H);
			break;
		case 3:
			printf("Please enter the element that you want to seach:");
			getchar();
			scanf("%d",&element);
			if(SearchHash(H,element,p,c)) printf("The element %d in %dth.\n",element,p);
			else printf("There is not the element.\n");
			break;
		case 4:
			ASL(H);
			break;
		case 5:
			printf("Please enter the element that you want to delete:");
			getchar();
			scanf("%d",&element);
			Delete(H,element);
			break;
		case 6:
			printf("Please enter the element that you want to insert:");
			getchar();
			scanf("%d",&element);
			InsertHash(H,element);
			break;
		case 0:
			exit(0);
		}
	}
}