#include"StackBracket.h"
#include"StackMaze.h"
#include"StackQueen.h"

void main(){
	int request,i,j;
	char c[100];
	Pos start,end;
	MazeType maze={{'#','#','#','#','#','#','#','#','#','#'},{'#','0','0','#','0','0','0','#','0','#'},{'#','0','0','#','0','0','0','#','0','#'},{'#','0','0','0','0','#','#','0','0','#'},{'#','0','#','#','#','0','0','0','0','#'},{'#','0','0','0','#','0','0','0','0','#'},{'#','0','#','0','0','0','#','0','0','#'},{'#','0','#','#','#','0','#','#','0','#'},{'#','#','0','0','0','0','0','0','0','#'},{'#','#','#','#','#','#','#','#','#','#'}};
	while(1){
		printf("*****Õ»µÄÓ¦ÓÃ*****\n");
		printf("0.Exit.\n");
		printf("1.Bracket application.\n");
		printf("2.Maze .\n");
		printf("3.The eight queens problem.\n");
		printf("Please input your request:");
		scanf("%d",&request);
		if(request==0) break;
		switch(request){
		case 1:
			printf("Please input your judgemental brackets:");
			getchar();
			gets(c);
			BracketsMatching(c);
			break;
		case 2:
			//printf("Please construct your maze by line(0 is pass;# is the wall.):\n");
			printf("Please input your entrance:");
			scanf("%d %d",&start.row,&start.column);
			printf("Please input your exit:");
			scanf("%d %d",&end.row,&end.column);
			//for(i=0;i<10;i++){
			//	printf("Please input a line:");
			//	for(j=0;j<10;j++) scanf("%c",&maze[i][j]);
			//}
			MazePath(maze,start,end);
			printf("The maze and path are shown below:\n");
			for(i=0;i<10;i++){
				for(j=0;j<10;j++){
					switch(*(maze[i]+j)){
					case '#':
						printf("#");
				        break;
					case 'X':
						printf(" ");
				        break;
				    case '0':
				        printf(" ");
				        break;
				    case '*':
				        printf("*");
				        break;
					}
				}
				printf("\n");
			}
			break;
		case 3:
			Queen();
			break;
		}
	}
}