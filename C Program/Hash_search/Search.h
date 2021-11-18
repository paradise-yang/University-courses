#include<stdio.h>
#include<stdlib.h>

//int hashsize[]={}
#define success 1
#define unsuccess 0
#define duplicate -1

typedef struct {
	int *elem;
	int count;
	int sizeindex;
}HashTable;

int InitDsTable(HashTable &H){
	int size,element=1,temp,i;
	printf("Please enter the length of the hashtable:");
	getchar();
	scanf("%d",&size);
	H.sizeindex=size;
	H.elem=(int *)malloc(size*sizeof(int));
	for(i=0;i<size;i++) H.elem[i]=0;
	printf("Pleas enter a elem:");
	getchar();
	scanf("%d",&element);
	while(element){
		//printf("Pleas enter a elem:");
		//scanf("%d",&element);
		temp=element%size;
		if(!H.elem[temp]) H.elem[temp]=element;
		else {
			while(H.elem[temp]) temp=(++temp)%size;
			H.elem[temp]=element;
		}
		H.count++;
		printf("Pleas enter a elem(Over is 0):");
		getchar();
		scanf("%d",&element);
	}
	return success;
}
int DestroyDSTable(HashTable &H){
	free(H.elem);
	return success;
}
int SearchHash(HashTable H,int key,int &p,int &c){
	p=key%H.sizeindex;
	c=1;
	while(H.elem[p] && key!=H.elem[p]){
		p=(++p)%H.sizeindex;
		c++;
	}
	if(key==H.elem[p]) return success;
	else return unsuccess;
}
int InsertHash(HashTable &H,int e){
	int c=0,p;
	if(SearchHash(H,e,p,c)) return duplicate;//´ý²åÈëÒÑ´æÔÚ
	else {
		H.elem[p]=e;
		++H.count;
		return success;
	}
}
int Traverse(HashTable H){
	int i;
	for(i=0;i<H.sizeindex;i++) printf("%-4d",i);
	printf("\n");
	for(i=0;i<H.sizeindex;i++){
		if(H.elem[i]) printf("%-4d",H.elem[i]);
		else printf("^   ");
	}
	printf("\n");
	return success;
}
int ASL(HashTable H){
	int i,p,c,su=0,fa=0;
	float success_ASL=0,fail_ASL=0;
	for(i=0;i<H.sizeindex;i++){
		if(H.elem[i]){
			SearchHash(H,H.elem[i],p,c);
			success_ASL+=c;
			++su;
		}
		else{
			SearchHash(H,H.elem[i],p,c);
			fail_ASL+=c;
			++fa;
		}
	}
	printf("The successful ASL=%.1f\n",success_ASL/su);
	printf("The unsuccessful ASL=%.1f\n",fail_ASL/fa);
	return success;
}
int Delete(HashTable &H,int key){
	int p,c;
	if(SearchHash(H,key,p,c)){
		H.elem[p]=0;
		H.count--;
	}
	else printf("There is not the element that you want to delete.\n");
	return success;
}
