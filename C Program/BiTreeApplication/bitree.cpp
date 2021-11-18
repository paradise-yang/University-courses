#include "BiTree.h"

void main(){
	char a[Max],*p[Max],tem,parent='#',lch='#',rch='#',brother='#';
	*p=a;
	BiTree *T=NULL,*copy=NULL,*insert=NULL;
	int request,nodes=0,leaves=0,depth,direction;
	while(1){
		printf("\n*****Ê÷µÄÓ¦ÓÃ*****\n");
		printf("0.Exit.\n");
		printf("1.Initialize the BiTree.\n");
		printf("2.The preorder sequence of the BiTree.\n");
		printf("3.The middle order sequence of the BiTree.\n");
		printf("4.The post sequence of the BiTree.\n");
		printf("5.The number of the nodes.\n");
		printf("6.The number of the leaves.\n");
		printf("7.The depth of the BiTree.\n");
		printf("8.The level traversal sequence of the BiTree.\n");
		printf("9.Search a concret node.\n");
		printf("10.Judge a BiTree.\n");
		printf("11.Exchange the left and right subtrees.\n");
		printf("12.Delete a subtree.\n");
		printf("13.Copy a BiTree.\n");
		printf("14.Insert a subtree.\n");
		printf("Please input your request:");
		scanf("%d",&request);
		if(!request) break;
		switch(request){
		case 1:
			printf("Please input a preorder extends a sequence:");
			getchar();
			gets(a);
			T=InilBiTree(T,p);
			break;
		case 2:
			PreOrder(T);
			printf("\n");
			break;
		case 3:
			InOrder(T);
			printf("\n");
			break;
		case 4:
			PostOrder(T);
			printf("\n");
			break;
		case 5:
			Nodes(T,nodes);
			printf("The number of the nodes of the BiTree is %d\n",nodes);
			break;
		case 6:
			Leaves(T,leaves);
			printf("The number of the leaves of the BiTree is %d\n",leaves);
			break;
		case 7:
			depth=Depth(T);
			printf("The depth of the BiTree is %d\n",depth);
			break;
		case 8:
			LayerTaversal(T);
			break;
		case 9:
			printf("Please input the node you want to seach:");
			getchar();
			scanf("%c",&tem);
			Search(T,tem,parent,lch,rch,brother);
			printf("The search results are as follows:\n");
			if(parent=='#') printf("The node doesn't have parent.");
			else {printf("The parent of the node is %c",parent);}
			printf("\n");
			if(lch=='#') printf("The node doesn't have a lchild.");
			else printf("The lchild of the node is %c",lch);
			printf("\n");
			if(rch=='#') printf("The node doesn't have a lchild.");
			else printf("The rchild of the node is %c",rch);
			printf("\n");
			if(brother=='#') printf("The node doesn't have a brother.");
			else printf("The brother of the node is %c",brother);
			printf("\n");
			printf("\n");
			break;
		case 10:
			if(Judge_BiSortTree(T)) printf("The BiTree is a binary sort tree.\n");
			else printf("The BiTree isn't a binary sort tree.\n");
			if(Judge_Complete_BiTree(T)) printf("The BiTree is a complete binary tree.\n");
			else printf("The BiTree isn't a complete binary tree.\n");
			if(Judge_AVL_BiTree(T)) printf("The BiTree is a balanced binary tree.\n");
			else printf("The BiTree isn't a balanced binary tree.\n");
			break;
		case 11:
			printf("Please input the node that you want to swap its subtrees:");
			getchar();
			scanf("%c",&tem);
			Exchange(T,tem);
			break;
		case 12:
			printf("Please enter which subtree of the node to delete(lchild is 1;rchild is 2):");
			getchar();
			scanf("%c %d",&tem,&direction);
			Delete(T,tem,direction);
			break;
		case 13:
			copy=CopyTree(T);
			printf("The copy is complted.\n");
			printf("The preorder sequence of duplicate version of the BiTree is ");
			PreOrder(copy);
			printf("\n");
			break;
		case 14:
			printf("Please enter the preorder extension sequence of the subtree to be inserted:");
			scanf("%s",a);
			*p=a;
			insert=InilBiTree(insert,p);
			printf("Please enter the which node that you want to insert:");
			getchar();
			scanf("%c",&tem);
			printf("Please enter which subtree of the node that you want to insert:");
			getchar();
			scanf("%d",&direction);
			Insert(T,insert,tem,direction);
			break;
		default :break;
		}
	}
}