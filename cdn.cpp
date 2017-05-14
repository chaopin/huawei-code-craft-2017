#include "deploy.h"
#include "lib_io.h"
#include "lib_time.h"
#include "stdio.h"

int main(int argc, char *argv[])
{
    //print_time("Begin");
    char *topo[MAX_EDGE_NUM];
    int line_num;
	argv[1]="./1 Intermediate/case6.txt";
	//argv[1]="./2 Advanced/case8.txt";//9有问题
	//argv[1]="./0 Primary/case1.txt";
	//argv[1]="testdouble.txt";
	//argv[1]="C:/Users/Administrator/Desktop/华为软件大赛/code/max_min/max_min/max_min_0323_linux/2diff/case1.txt";
	//argv[2]="output.txt";
	argv[2]="G3path.txt";
    char *topo_file = argv[1];

    line_num = read_file(topo, MAX_EDGE_NUM, topo_file);

    printf("line num is :%d \n", line_num);
    if (line_num == 0)
    {
        printf("Please input valid topo file.\n");
        return -1;
    }

    char *result_file = argv[2];

    deploy_server(topo, line_num, result_file);

    release_buff(topo, line_num);

   // print_time("End");

	return 0;
}

