#include"Graph.h"

void main(){
	Graph G,Gn;
	int request;
	char v,v1,v2;
	printf("\n*****图的应用*****\n");
	printf("0.Exit.\n");
	printf("1.Initialize the Graph.\n");
	printf("2.DFS.\n");
	printf("3.Insert a Arc.\n");
	printf("4.Delete a Arc.\n");
	printf("5.Calculate the degree of a vertex.\n");
	printf("6.BFS.\n");
	printf("7.The shortest path pf two vertexs.\n");
	printf("8.Minimum cost spanning tree(Prim).\n");
	while(1){
		//printf("\n*****图的应用*****\n");
		//printf("0.Exit.\n");
		//printf("1.Initialize the Graph.\n");
		//printf("2.DFS.\n");
		printf("\nPlease input your request:");
		scanf("%d",&request);
		if(!request) break;
		switch(request){
		case 1:
			CreateUDN(G);
			break;
		case 2:
			DFSTraverse(G);
			break;
		case 3:
			InsertArc(G);
			break;
		case 4:
			DeleteArc(G);
			break;
		case 5:
			printf("Please enter the vertex that you want to calculate its degree:");
			getchar();
			scanf("%c",&v);
			DegreeVex(G,v);
			break;
		case 6:
			BFSTraverse(G);
			break;
		case 10:
			printf("Please enter the two vertexs:");
			getchar();
			scanf("%c%c",&v1,&v2);
			MinPath(G,v1,v2);
			break;
		case 7:
			CreateUDN_1(Gn);
			printf("Please enter the root vertex:");
			getchar();
			scanf("%c",&v);
			Dijkstra(Gn,v);
			break;
		case 8:
			CreateUDN_1(Gn);
			printf("Please enter the root vertex:");
			getchar();
			scanf("%c",&v);
			Prim(Gn,v);
			break;
		//case 9:
		//	ShortestPath(Gn);
		//	break;
		default: break;
		}
	}
}