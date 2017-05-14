
#include "deploy.h"
#include  <stdio.h>
#include <malloc.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <queue>
#include <ctime>
#include <cmath>


#include <time.h>
#include <memory.h>
#include <stack>
#include <string>
#include <set>
//#define FILENODE 0
using namespace std;



#define MAX_EDGE 10000
#define MAX_VER 1000
#define MAX_DIS  50000
#define MAX_CST  50000
#define MAX_COS  50000

#define M 65535
#define RELOAD 1
#define MERGE 1
#define DELETE 1
#define SIMPLEXVER 1
int cntResi=0;
//#define COUTINFO	1	
typedef struct
{
	unsigned short ind_s, ind_e;
	unsigned short limit;
	float cost;
	//the first two guys are current value,the middle two are backupvalue,the last two are best
	float cost_dsg[6];
	unsigned short resi_dsg[6];
	unsigned short bnd_dsg[6];
	unsigned short path;
	bool flag_chg;
}edge, *pe;
typedef struct
{
	unsigned short *edge;
	unsigned short *oppo;
	bool *edge_flag;

	unsigned short base_cost;
	unsigned short cons_id;
	unsigned short num_edge;
	bool server, customer;
	unsigned short demand;
	short given;
	float per_cus;
	float judge;
	float judge_distance;
	unsigned short num_visited;

	float dist;
	bool flag_q;
	unsigned short pre;
	unsigned short pre_arc;
	unsigned short total_output;
	short ser_level;
}vertex, *pv;
typedef struct pri_que
{
	unsigned short id;
	unsigned short demands;
	float percost;
	float cost;
	friend bool operator<(pri_que p1, pri_que p2)
	{
		if (p1.percost>p2.percost)
			return p1.percost>p2.percost;
		else if (p1.percost == p2.percost)
		{
			return p1.demands<p2.demands;
		}
		else
			return p1.percost>p2.percost;
	}
}pri_que;
typedef struct
{
	unsigned short ser_upper[11];
	unsigned short ser_cost[11];
	unsigned short ser_num;

	double dist;
	short pre;
	short band;
	unsigned int total_consume;

	unsigned short *server;
	float *distance;
	unsigned short *list_ser;
	unsigned int *flow;
	unsigned int *flow_bkup;
	int *given_bkup;
	unsigned int *flow_best;
	int *given_best;
	short *serLevbkup;
	short *serLevbest;

	unsigned int opt;
	short *solveNode;
	unsigned short *solveLink;
	unsigned short *solveSerRemain;

	short *opt_solve;
	unsigned int *opt_flow;
	unsigned short *opt_level;
	short *flag_black;
}se, *pse;
typedef struct
{
	unsigned short top, end, length, cntn;
	unsigned short *queue;
}que, *que_p;
struct Problem
{
	unsigned short numLink;
	unsigned short numNetNode;
	unsigned short numConNode;
	unsigned short numServerLevel;		// Number of server's levels
	unsigned short serverAbility[11];		// The ability of server
	unsigned short serverCost[11];			// The cost of server
	unsigned short *nodeCost;				// The cost of building a server at the node
	unsigned short *linkIndex1;
	unsigned short *linkIndex2;
	unsigned short *linkBandwidth;
	unsigned short *linkCost;
	short *nodeNeed;
};
struct Change
{
	unsigned short changeNum;
	unsigned short changeIndex[100];
	short changeOld[100];
	short changeNew[100];
};


//unsigned short num_nn,num_link,num_cn,cost_ser;
unsigned short num_nn, num_link, num_cn;
que_p min_max_q;
unsigned short temp_pop, ind, temp_visit, ind_e;
unsigned short v_s, v_e;
unsigned char temp_path;
float temp_cost;
unsigned int consume;
short temp_ind, temp_arc;
unsigned char temp_p;
Problem prob;
unsigned long potential_cost;
float alpha;
clock_t startTime, endTime;


inline bool check_opt(pse sta_end, unsigned int opt_potential)
{
	if (opt_potential < sta_end->opt)
	{
		sta_end->opt = opt_potential;
		for (int i = 0; i<num_link * 2; i++) sta_end->opt_solve[i] = sta_end->solveLink[i];
		for (int i = num_link * 2; i<num_nn + num_link * 2; i++) sta_end->opt_solve[i] = sta_end->solveNode[i - num_link * 2];
		return true;
	}
	return false;
}
inline bool jump_time(unsigned int time_out)
{
	endTime = clock();
	if (((double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC) > time_out) return true;
	return false;
}
inline void Black(const Problem *prob, short *black);
bool nonempty_q(que_p queue);
que_p init_q(unsigned short length);
inline void initial(Problem *prob);
inline void restore(pv list_v, pe list_e, pse sta_end);

inline unsigned int simplex(const Problem *prob, const short *solveNode, unsigned short *solveLink, unsigned short *serverRemaining);
inline unsigned int simplex_simple(const Problem *prob, const short *solveNode);
inline int min_max(pv list_v, pe list_e, pse stt_end, int total, unsigned long input, short *cost);

inline int road_milp(pv list_v, pe list_e, pse sta_end, unsigned long totalcoms, char *output, short *solve);
inline int cal_cost(pse sta_end, unsigned short ser_flow);
inline int cal_level(pse sta_end, unsigned short ser_flow);
inline void adj_ansMid(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out);
inline void adj_ansLge(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out);
inline int iter_ansMid(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out);
inline int iter_ansLge(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out);

inline void adj_ansMid2(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out);
inline void adj_ansLge2(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out);
inline int iter_ansMid2(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out);
inline int iter_ansLge2(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out);
inline void DelRelAdj(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned int time_out);
inline void cnvtGraph(pv list_v, pe list_e, pse sta_end);
inline unsigned int resi_map(pv list_v, pe list_e, pse sta_end, unsigned short *solveLink, const Change *change);
inline void backup_graph(pv list_v, pe list_e, pse sta_end);
inline void recover_graph(pv list_v, pe list_e, pse sta_end);
inline void bkup_best(pv list_v, pe list_e, pse sta_end);
inline void recv_best(pv list_v, pe list_e, pse sta_end);

void deploy_server(char * topo[MAX_EDGE_NUM], int line_num, char * filename)
{
	unsigned long totalcoms = 0;

	startTime = clock();
	sscanf(topo[0], "%hd %hd %hd\n", &num_nn, &num_link, &num_cn);

	prob.numLink = num_link;	prob.numNetNode = num_nn;	prob.numConNode = num_cn;
	prob.linkIndex1 = (unsigned short*)malloc(num_link*sizeof(unsigned short));
	prob.linkIndex2 = (unsigned short*)malloc(num_link*sizeof(unsigned short));
	prob.linkBandwidth = (unsigned short*)malloc(num_link*sizeof(unsigned short));
	prob.linkCost = (unsigned short*)malloc(num_link*sizeof(unsigned short));
	prob.nodeNeed = (short*)malloc(num_nn*sizeof(short));

	pe list_e = (pe)malloc((num_link)*sizeof(edge));	pv list_v = (pv)malloc((num_nn)*sizeof(vertex));
	pse sta_end = (pse)malloc(sizeof(se));
	memset(sta_end, 0, sizeof(se));
	memset(list_e, 0, (num_link)*sizeof(edge));
	memset(list_v, 0, (num_nn)*sizeof(vertex));

	sta_end->server = (unsigned short *)malloc((num_nn)*sizeof(unsigned short));
	sta_end->list_ser = (unsigned short *)malloc((num_nn)*sizeof(unsigned short));
	sta_end->flow = (unsigned int*)malloc((num_nn)*sizeof(unsigned int));
	sta_end->flow_bkup = (unsigned int*)malloc((num_nn)*sizeof(unsigned int));
	sta_end->given_bkup = (int*)malloc((num_nn)*sizeof(int));
	sta_end->flow_best = (unsigned int*)malloc((num_nn)*sizeof(unsigned int));
	sta_end->given_best = (int*)malloc((num_nn)*sizeof(int));
	sta_end->serLevbkup = (short *)malloc((num_nn)*sizeof(short));
	sta_end->serLevbest = (short *)malloc((num_nn)*sizeof(short));
	sta_end->distance = (float *)malloc((num_nn)*sizeof(float));
	sta_end->opt_solve = (short *)malloc((num_nn + 2 * num_link)*sizeof(short));
	sta_end->opt_flow = (unsigned int *)malloc((num_nn)*sizeof(unsigned int));
	sta_end->opt_level = (unsigned short *)malloc((num_nn)*sizeof(unsigned short));

	sta_end->flag_black = (short *)malloc((num_nn)*sizeof(short));
	sta_end->solveLink = (unsigned short*)malloc((num_link * 2)*sizeof(unsigned short));
	sta_end->solveNode = (short*)malloc(num_nn*sizeof(short));
	sta_end->solveSerRemain = (unsigned short*)malloc((num_nn)*sizeof(unsigned short));

	min_max_q = init_q(num_nn + 2);

	memset(sta_end->solveNode, -1, sizeof(short)* num_nn);

	unsigned short ser_level, ser_upper, ser_cost;
	unsigned short temp_row = 2;
	sta_end->ser_num = 0;
	int temp_ser_num = 0;
	while (1)
	{
		ser_level = 0; ser_upper = 0; ser_cost = 0;
		sscanf(topo[temp_row], "%hd %hd %hd\n", &ser_level, &ser_upper, &ser_cost);
		prob.serverAbility[sta_end->ser_num] = sta_end->ser_upper[sta_end->ser_num] = ser_upper;
		prob.serverCost[sta_end->ser_num] = sta_end->ser_cost[sta_end->ser_num] = ser_cost;
		if (ser_level == 0 && ser_upper == 0 && ser_cost == 0) break;
		temp_row++;
		sta_end->ser_num++;
	}
	prob.numServerLevel = sta_end->ser_num;
	unsigned short ind_start, ind_end, temp_c;
	prob.nodeCost = (unsigned short*)malloc(num_nn*sizeof(unsigned short));
	for (unsigned short i = temp_row + 1; i<temp_row + 1 + num_nn; i++)
	{
		sscanf(topo[i], "%hd %hd\n", &ser_level, &list_v[i - temp_row - 1].base_cost);
		prob.nodeCost[i - temp_row - 1] = list_v[i - temp_row - 1].base_cost;
	}
	//cal the num of one vertex
	unsigned short index;
	for (ind = 0; ind<num_link; ind++)
	{
		index = temp_row + 2 + num_nn + ind;
		sscanf(topo[index], "%hd %hd %hd %hd\n", &list_e[ind].ind_s, &list_e[ind].ind_e, &list_e[ind].limit, &temp_c);
		prob.linkIndex1[ind] = list_e[ind].ind_s;
		prob.linkIndex2[ind] = list_e[ind].ind_e;
		prob.linkBandwidth[ind] = list_e[ind].limit;
		prob.linkCost[ind] = temp_c;

		list_e[ind].cost = temp_c;
		list_e[ind].cost_dsg[0] = list_e[ind].cost_dsg[1] = list_e[ind].cost;
		list_e[ind].resi_dsg[0] = list_e[ind].resi_dsg[1] = list_e[ind].limit;
		list_e[ind].bnd_dsg[0] = list_e[ind].bnd_dsg[1] = list_e[ind].limit;

		ind_start = list_e[ind].ind_s;
		ind_end = list_e[ind].ind_e;

		list_v[ind_start].num_edge++;
		list_v[ind_end].num_edge++;
	}
	//reset the graph
	for (ind = 0; ind<num_nn; ind++)
	{
		list_v[ind].edge = (unsigned short *)malloc(list_v[ind].num_edge*sizeof(unsigned short));
		list_v[ind].oppo = (unsigned short *)malloc(list_v[ind].num_edge*sizeof(unsigned short));
		list_v[ind].edge_flag = (bool*)malloc(list_v[ind].num_edge*sizeof(bool));
		list_v[ind].num_edge = 0;
	}
	//read in the edge
	for (ind = 0; ind<num_link; ind++)
	{
		ind_start = list_e[ind].ind_s;
		ind_end = list_e[ind].ind_e;

		list_v[ind_start].edge[list_v[ind_start].num_edge] = ind;
		list_v[ind_start].oppo[list_v[ind_start].num_edge] = ind_end;
		list_v[ind_start].edge_flag[list_v[ind_start].num_edge] = 0;
		list_v[ind_start].num_edge++;
		list_v[ind_start].total_output += list_e[ind].limit;

		list_v[ind_end].edge[list_v[ind_end].num_edge] = ind;
		list_v[ind_end].oppo[list_v[ind_end].num_edge] = ind_start;
		list_v[ind_end].edge_flag[list_v[ind_end].num_edge] = 1;
		list_v[ind_end].num_edge++;
		list_v[ind_end].total_output += list_e[ind].limit;
	}
	//read consume
	unsigned short id_cus, id_net, demand;
	unsigned short start_consume;
	unsigned short sercost;
	memset(prob.nodeNeed, 0, (num_nn)*sizeof(short));
	for (ind = 0; ind<num_cn; ind++)
	{
		start_consume = temp_row + 3 + num_nn + num_link + ind;
		sscanf(topo[start_consume], "%hd %hd %hd\n", &id_cus, &id_net, &demand);
		prob.nodeNeed[id_net] = -demand;
		totalcoms += demand;
		list_v[id_net].cons_id = id_cus;
		list_v[id_net].demand = demand;
		list_v[id_net].customer = true;

		sercost = cal_cost(sta_end, demand);
		list_v[id_net].per_cus = (float)(sercost + list_v[id_net].base_cost) / (float)demand;
		list_v[id_net].num_visited++;
		list_v[id_net].total_output += list_v[id_net].demand;
	}
	sta_end->total_consume = totalcoms;
	sta_end->opt = (2000 + sta_end->ser_cost[sta_end->ser_num - 1])*num_cn;
	Black(&prob, sta_end->flag_black);

	unsigned int time_out = 85000;
	unsigned int time_upper = 89000;
	unsigned int temp_opt;
	int MidSele = 1;
	int LgeSele = 1;
	if (num_nn<700)
	{
		if (MidSele == 1)
		{
			time_out = 15000;
			iter_ansMid(list_v, list_e, sta_end, &prob, totalcoms, time_out);
			time_out = time_upper;
			cnvtGraph(list_v, list_e, sta_end);
			while (1)
			{
				while (1)
				{
					temp_opt = sta_end->opt;
					adj_ansMid(list_v, list_e, sta_end, &prob, totalcoms, time_out);
					endTime = clock(); if (((double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC)>time_out) goto whole_end;
					if (temp_opt == sta_end->opt) break;
				}
				endTime = clock(); if (((double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC)>time_out) goto whole_end;
				DelRelAdj(list_v, list_e, sta_end, &prob, time_out);
				endTime = clock(); if (((double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC)>time_out) goto whole_end;
			}
		}
		else
		{
			time_out = 15000;
			iter_ansMid2(list_v, list_e, sta_end, &prob, totalcoms, time_out);
			time_out = time_upper;
			cnvtGraph(list_v, list_e, sta_end);
			while (1)
			{
				while (1)
				{
					temp_opt = sta_end->opt;
					adj_ansMid2(list_v, list_e, sta_end, &prob, totalcoms, time_out);
					endTime = clock(); if (((double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC)>time_out) goto whole_end;
					if (temp_opt == sta_end->opt) break;
				}
				endTime = clock(); if (((double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC)>time_out) goto whole_end;
				DelRelAdj(list_v, list_e, sta_end, &prob, time_out);
				endTime = clock(); if (((double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC)>time_out) goto whole_end;
			}

		}
	}
	else
	{
		if (LgeSele == 1)
		{
			time_out = 15000;
			iter_ansLge(list_v, list_e, sta_end, &prob, totalcoms, time_out);
			time_out = time_upper;
			cnvtGraph(list_v, list_e, sta_end);
			while (1)
			{
				adj_ansLge(list_v, list_e, sta_end, &prob, totalcoms, time_out);
				endTime = clock(); if (((double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC)>time_out) goto whole_end;
			}
		}
		else
		{
			time_out = 15000;
			iter_ansLge2(list_v, list_e, sta_end, &prob, totalcoms, time_out);
			time_out = time_upper;
			cnvtGraph(list_v, list_e, sta_end);
			while (1)
			{
				adj_ansLge2(list_v, list_e, sta_end, &prob, totalcoms, time_out);
				endTime = clock(); if (((double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC)>time_out) goto whole_end;
			}
		}
	}
whole_end:

	endTime = clock();

#ifdef SIMPLEXVER
	recv_best(list_v, list_e, sta_end);
	for (int i = 0; i<num_nn; i++)
	{
		sta_end->solveNode[i] = -1;
		if (sta_end->flow[i]>0)
		{
			sta_end->solveNode[i] = cal_level(sta_end, sta_end->flow[i]);
		}
	}
	cout << "simplex is " << simplex(&prob, sta_end->solveNode, sta_end->solveLink, sta_end->solveSerRemain);
	cout << "opt is " << sta_end->opt << endl;
	cout<<"cntResi is "<<cntResi<<endl;
	cout<<((double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC)<<" ms"<<endl;
#endif
	char *topo_file = (char*)malloc(200000 * sizeof(char));
	restore(list_v, list_e, sta_end);
	road_milp(list_v, list_e, sta_end, totalcoms, topo_file, sta_end->opt_solve);
	write_result(topo_file, filename);
	free(topo_file);

	free(list_e);
	free(list_v);

	free(sta_end->server);
	free(sta_end->list_ser);
	free(sta_end->flow);
	free(sta_end->distance);
	free(sta_end->opt_solve);
	free(sta_end->opt_flow);
	free(sta_end->opt_level);
	free(sta_end->solveLink);
	free(sta_end->solveNode);
	free(sta_end->solveSerRemain);
	free(sta_end);
}
void Black(const Problem *prob, short *black)
{
	unsigned short *degree = new unsigned short[prob->numNetNode];
	memset(degree, 0, sizeof(short)* prob->numNetNode);
	unsigned short *outCapability = new unsigned short[prob->numNetNode];
	memset(outCapability, 0, sizeof(short)* prob->numNetNode);
	for (int i = 0; i < prob->numLink; i++)
	{
		degree[prob->linkIndex1[i]]++;
		degree[prob->linkIndex2[i]]++;
		outCapability[prob->linkIndex1[i]] += prob->linkBandwidth[i];
		outCapability[prob->linkIndex2[i]] += prob->linkBandwidth[i];
	}
	unsigned short *outAndneed = new unsigned short[prob->numNetNode];
	memset(outAndneed, 0, sizeof(short)* prob->numNetNode);
	for (int i = 0; i < prob->numNetNode; i++)
		outAndneed[i] = (short)outCapability[i] - prob->nodeNeed[i];
	unsigned short *bestServer = new unsigned short[prob->numNetNode];
	memset(bestServer, 0, sizeof(short)* prob->numNetNode);
	for (int i = 0; i < prob->numNetNode; i++)
	{

		unsigned short level = 0;
		double bestServerCost = 65535;
		double thisServerCost = 0;
		while (outAndneed[i] > prob->serverAbility[level])
		{
			thisServerCost = (double)(prob->serverCost[level] + prob->nodeCost[i]) / prob->serverAbility[level];
			if (thisServerCost < bestServerCost)
			{
				bestServerCost = thisServerCost;
				bestServer[i] = level;
			}
			level++;
			if (level == prob->numServerLevel)break;
		}
		if (level != prob->numServerLevel)
		{
			thisServerCost = (double)(prob->serverCost[level] + prob->nodeCost[i]) / outAndneed[i];
			if (thisServerCost < bestServerCost)
			{
				bestServerCost = thisServerCost;
				bestServer[i] = level;
			}
		}
	}
	for (int i = 0; i < prob->numNetNode; i++)
	{
		if (degree[i] != 1)
			black[i] = bestServer[i];
		else
			black[i] = -1;
	}
	delete[] degree;
	delete[] outCapability;
	delete[] outAndneed;
	delete[] bestServer;
}
bool nonempty_q(que_p queue)
{
	/*if(queue->cntn) return true;
	else return false;*/
	return (bool)queue->cntn;
}
que_p init_q(unsigned short length)
{
	que_p queue = (que_p)malloc(sizeof(que));
	queue->top = 0;
	queue->end = 0;
	queue->cntn = 0;
	queue->length = length;
	queue->queue = (unsigned short*)malloc(length*sizeof(unsigned short));
	return queue;
}
void initial(Problem *prob)
{
	ifstream infile("case1.txt");
	string tempString;
	getline(infile, tempString, ' ');		prob->numNetNode = stoi(tempString);
	getline(infile, tempString, ' ');		prob->numLink = stoi(tempString);
	getline(infile, tempString);			prob->numConNode = stoi(tempString);
	getline(infile, tempString);
	/* Read the Levels of Server */
	prob->numServerLevel = 6;
	for (int i = 0; i < prob->numServerLevel; i++)
	{
		getline(infile, tempString, ' ');
		getline(infile, tempString, ' ');
		prob->serverAbility[i] = stoi(tempString);
		getline(infile, tempString);
		prob->serverCost[i] = stoi(tempString);
	}
	getline(infile, tempString);
	/* Read the Cost of Building a Server at the Node */
	prob->nodeCost = new unsigned short[prob->numNetNode];
	for (int i = 0; i < prob->numNetNode; i++)
	{
		getline(infile, tempString, ' ');
		getline(infile, tempString);
		prob->nodeCost[i] = stoi(tempString);
	}
	getline(infile, tempString);
	/* Definition of Variables of Link */
	prob->linkIndex1 = new unsigned short[prob->numLink];
	prob->linkIndex2 = new unsigned short[prob->numLink];
	prob->linkBandwidth = new unsigned short[prob->numLink];
	prob->linkCost = new unsigned short[prob->numLink];
	/* Read Links */
	for (int i = 0; i < prob->numLink; i++)
	{
		getline(infile, tempString, ' ');		prob->linkIndex1[i] = stoi(tempString);
		getline(infile, tempString, ' ');		prob->linkIndex2[i] = stoi(tempString);
		getline(infile, tempString, ' ');		prob->linkBandwidth[i] = stoi(tempString);
		getline(infile, tempString);			prob->linkCost[i] = stoi(tempString);
	}
	getline(infile, tempString);
	/* Read Nodes */
	prob->nodeNeed = new short[prob->numNetNode];
	memset(prob->nodeNeed, 0, sizeof(short)* prob->numNetNode);
	unsigned short nodeNet;	//	The Network Node that the Consumer Node connected
	unsigned short conNeed;	//	The Network Node that the Consumer Node connected
	for (int i = 0; i < prob->numConNode; i++)
	{
		getline(infile, tempString, ' ');
		getline(infile, tempString, ' ');		nodeNet = stoi(tempString);
		getline(infile, tempString);			conNeed = stoi(tempString);
		prob->nodeNeed[nodeNet] = 0 - conNeed;
	}
	//infile.close();
}
void restore(pv list_v, pe list_e, pse sta_end)
{
	//back to the original
	//clear the edge
	for (ind = 0; ind<num_link; ind++)
	{
		//if(!list_e[ind].flag_chg) continue;
		list_e[ind].flag_chg = false;
		list_e[ind].cost_dsg[0] = list_e[ind].cost_dsg[1] = list_e[ind].cost;
		list_e[ind].resi_dsg[0] = list_e[ind].resi_dsg[1] = list_e[ind].limit;
		list_e[ind].bnd_dsg[0] = list_e[ind].bnd_dsg[1] = list_e[ind].limit;
	}
	//clear the start pnt
	for (ind = 0; ind<num_nn; ind++)
	{
		sta_end->distance[ind] = MAX_DIS;
		sta_end->server[ind] = 0;
		sta_end->flow[ind] = 0;
		if (!list_v[ind].customer) continue;
		list_v[ind].given = 0;
	}
}

int iter_ansMid(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out)
{
	unsigned short ind;
	unsigned int temp_cost_sum, temp_need_sum;
	short temp_cost;
	priority_queue<pri_que> p_quque;
	pri_que temp_pri;
	unsigned long  temp_totalcoms = totalcoms;
	stack<unsigned short> stack_delete;
	unsigned int cost_simp;
	float alpha = (float)1.2;
	float alpha_s, alpha_e, alpha_step;
	int temp_ser_level;
	alpha_s = (float)0.5; alpha_e = (float)1.6; alpha_step = (float)0.4;

	alpha = alpha_s;

	for (ind = 0; ind<num_nn; ind++)
	{
		list_v[ind].ser_level = sta_end->flag_black[ind];
		if (sta_end->flag_black[ind] >= 0)
		{
			temp_cost_sum = list_v[ind].base_cost;
			temp_need_sum = list_v[ind].demand;
			temp_ser_level = sta_end->flag_black[ind];
			for (int edge_ind = 0; edge_ind<list_v[ind].num_edge; edge_ind++)
			{
				int edge = list_v[ind].edge[edge_ind];
				temp_cost_sum += (unsigned int)list_e[edge].cost*list_e[edge].limit;
				temp_need_sum += list_e[edge].limit;
			}
			list_v[ind].judge = float(temp_cost_sum + sta_end->ser_cost[temp_ser_level]) / float(temp_need_sum);
		}
		else list_v[ind].judge = MAX_DIS;
	}

	ind = 0;
	while (alpha<alpha_e)
	{
		for (ind = 0; ind<num_link; ind++)
		{
			list_e[ind].cost_dsg[0] = list_e[ind].cost_dsg[1] = list_e[ind].cost;
			list_e[ind].resi_dsg[0] = list_e[ind].resi_dsg[1] = list_e[ind].limit;
			list_e[ind].bnd_dsg[0] = list_e[ind].bnd_dsg[1] = list_e[ind].limit;
		}
		for (ind = 0; ind<num_nn; ind++)
		{
			sta_end->flow[ind] = 0;
			list_v[ind].given = 0;
			if (list_v[ind].ser_level >= 0)
			{
				sta_end->server[ind] = true;
				sta_end->distance[ind] = alpha*list_v[ind].judge;
			}
			else
			{
				sta_end->server[ind] = false;
				sta_end->distance[ind] = MAX_DIS;
			}
		}
		if (jump_time(time_out)) goto whole_end;
		totalcoms = temp_totalcoms;
		min_max(list_v, list_e, sta_end, num_nn, totalcoms, &temp_cost);
		for (int i = 0; i<num_nn; i++)
		{
			sta_end->solveNode[i] = -1;//solveNode[i]=0;
			if (sta_end->flow[i]>0) sta_end->solveNode[i] = sta_end->flag_black[i];
		}
		if (jump_time(time_out)) goto whole_end;
		cost_simp = simplex(prob, sta_end->solveNode, sta_end->solveLink, sta_end->solveSerRemain);
		check_opt(sta_end, cost_simp);
		alpha = alpha + alpha_step;
		if (jump_time(time_out)) goto whole_end;
	}

whole_end:
	return 0;
}
void adj_ansMid(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out)
{

	priority_queue<pri_que> chg_queue;
	pri_que temp_pri;
	Change change;
	unsigned int sum_cst, sum_need, opt;
	short temp_top, oppo_temp, temp_ser_level;
	unsigned short oppo_temp_level, flowOppo, flowTop;

	recv_best(list_v, list_e, sta_end);

	for (ind = 0; ind<num_nn; ind++)
	{
		if (sta_end->flow[ind]>0)
		{
			sum_cst = list_v[ind].base_cost + sta_end->ser_cost[list_v[ind].ser_level];
			sum_need = sta_end->flow[ind];
			temp_pri.id = ind;
			temp_pri.demands = list_v[ind].demand;
			temp_pri.percost = -float(sum_cst) / float(sum_need);
			temp_pri.cost = -float(sum_cst) / float(sum_need);
			chg_queue.push(temp_pri);
		}
	}
	if (jump_time(time_out)) { goto whole_end; }
	while (!chg_queue.empty())
	{
		recv_best(list_v, list_e, sta_end);
		temp_pri = chg_queue.top();
		chg_queue.pop();
		temp_top = temp_pri.id;
		temp_ser_level = list_v[temp_top].ser_level;
		//delete current solution
		change.changeNum = 1;
		change.changeIndex[0] = temp_top;
		change.changeOld[0] = temp_ser_level;
		change.changeNew[0] = -1;
		opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
		if (check_opt(sta_end, opt)) {
#ifdef COUTINFO				
			cout << "Delete simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
			cout << "Delete opt" << opt << endl;
#endif
			bkup_best(list_v, list_e, sta_end);
		}
		change.changeIndex[1] = temp_top;
		change.changeOld[1] = temp_ser_level;
		change.changeNew[1] = -1;
		change.changeNum = 2;
		for (int chg_ind = 0; chg_ind<list_v[temp_top].num_edge; chg_ind++)
		{
			recover_graph(list_v, list_e, sta_end);
			oppo_temp = list_v[temp_top].oppo[chg_ind];
			change.changeIndex[0] = oppo_temp;
			change.changeOld[0] = list_v[oppo_temp].ser_level;
			change.changeNew[1] = -1;
			if (list_v[oppo_temp].ser_level<0 && sta_end->flag_black[oppo_temp] >= 0)
			{
				change.changeNew[0] = temp_ser_level;
#ifdef COUTINFO				
				cout << "move start" << endl;
#endif
				opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
				if (check_opt(sta_end, opt))
				{
#ifdef COUTINFO				
					cout << "simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
					cout << "Move opt" << opt << endl;
#endif
					bkup_best(list_v, list_e, sta_end);
				}
				if (jump_time(time_out)) goto whole_end;
			}
			else if (list_v[oppo_temp].ser_level >= 0 && sta_end->flag_black[oppo_temp] >= 0)
			{
#ifdef MERGE
				oppo_temp_level = list_v[oppo_temp].ser_level;
				flowOppo = sta_end->flow[oppo_temp];
				flowTop = sta_end->flow[temp_top];
				//cout<<"Merge start"<<endl;
				if ((flowOppo + flowTop)>sta_end->ser_upper[sta_end->flag_black[oppo_temp]])
				{
					change.changeNew[0] = sta_end->flag_black[oppo_temp];
					change.changeNew[1] = cal_level(sta_end, (flowOppo + flowTop) - sta_end->ser_upper[sta_end->flag_black[oppo_temp]]);
				}
				else
				{
					change.changeNew[0] = cal_level(sta_end, (flowOppo + flowTop));
					change.changeNew[1] = -1;
				}
				list_v[oppo_temp].ser_level = change.changeNew[0];
				list_v[temp_top].ser_level = change.changeNew[1];
				opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
				if (check_opt(sta_end, opt))
				{
#ifdef COUTINFO				
					cout << "simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
					cout << "MERGE opt" << opt << endl;
#endif
					bkup_best(list_v, list_e, sta_end);
				}
				if (jump_time(time_out))  goto whole_end;
#endif         
			}
		}
		recv_best(list_v, list_e, sta_end);
#ifdef RELOAD
		short idMaxresi, idMinresi, resiMax, resiMin, resiTemp;
		resiMax = 0; resiMin = sta_end->ser_upper[sta_end->ser_num - 1];
		idMaxresi = -1; idMinresi = -1;
		for (int i = 0; i<num_nn; i++)
		{
			if (sta_end->flow[i])
			{
				resiTemp = sta_end->ser_upper[list_v[i].ser_level] - sta_end->flow[i];
				if (resiTemp<resiMin && (list_v[i].ser_level<sta_end->ser_num - 1))
					//if(resiTemp<resiMin&&list_v[i].ser_level<sta_end->flag_black[i]) 
				{
					idMinresi = i;
					resiMin = resiTemp;
				}
				if (resiTemp>resiMax)
				{
					idMaxresi = i;
					resiMax = resiTemp;
				}
			}
		}
		if (idMaxresi != -1 && idMinresi != -1 && idMinresi != idMaxresi)
		{
			change.changeNum = 2;
			change.changeIndex[0] = idMaxresi;
			change.changeOld[0] = list_v[idMaxresi].ser_level;
			change.changeNew[0] = list_v[idMaxresi].ser_level - 1;
			change.changeIndex[1] = idMinresi;
			change.changeOld[1] = list_v[idMinresi].ser_level;
			change.changeNew[1] = list_v[idMinresi].ser_level + 1;
#ifdef COUTINFO				
			cout << "Reload start" << endl;
#endif
			opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
			if (check_opt(sta_end, opt))
			{
#ifdef COUTINFO				
				cout << "simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
				cout << "RELOAD opt" << opt << endl;
#endif
				bkup_best(list_v, list_e, sta_end);
			}
			if (jump_time(time_out)) goto whole_end;
		}
		recv_best(list_v, list_e, sta_end);
		if (jump_time(time_out)) goto whole_end;
#endif
		backup_graph(list_v, list_e, sta_end);
	}
whole_end:;
}
int iter_ansLge(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out)
{
	unsigned short ind;
	unsigned int temp_cost_sum, temp_need_sum;
	short temp_cost;
	priority_queue<pri_que> p_quque;
	pri_que temp_pri;
	unsigned long  temp_totalcoms = totalcoms;
	stack<unsigned short> stack_delete;
	unsigned int cost_simp;
	float alpha = (float)1.2;
	float alpha_s, alpha_e, alpha_step;
	int temp_ser_level;
	alpha_s = (float)0.5; alpha_e = (float)1.6; alpha_step = (float)0.2;

	alpha = alpha_s;

	for (ind = 0; ind<num_nn; ind++)
	{
		list_v[ind].ser_level = sta_end->flag_black[ind];
		if (sta_end->flag_black[ind] >= 0)
		{
			temp_cost_sum = list_v[ind].base_cost;
			temp_need_sum = list_v[ind].demand;
			temp_ser_level = sta_end->flag_black[ind];
			for (int edge_ind = 0; edge_ind<list_v[ind].num_edge; edge_ind++)
			{
				int edge = list_v[ind].edge[edge_ind];
				temp_cost_sum += (unsigned int)list_e[edge].cost*list_e[edge].limit;
				temp_need_sum += list_e[edge].limit;
			}
			list_v[ind].judge = float(temp_cost_sum + sta_end->ser_cost[temp_ser_level]) / float(temp_need_sum);
		}
		else list_v[ind].judge = MAX_DIS;
	}

	ind = 0;
	for (ind = 0; ind<num_link; ind++)
	{
		list_e[ind].cost_dsg[0] = list_e[ind].cost_dsg[1] = list_e[ind].cost;
		list_e[ind].resi_dsg[0] = list_e[ind].resi_dsg[1] = list_e[ind].limit;
		list_e[ind].bnd_dsg[0] = list_e[ind].bnd_dsg[1] = list_e[ind].limit;
	}
	for (ind = 0; ind<num_nn; ind++)
	{
		sta_end->flow[ind] = 0;
		list_v[ind].given = 0;
		if (list_v[ind].ser_level >= 0)
		{
			sta_end->server[ind] = true;
			sta_end->distance[ind] = alpha*list_v[ind].judge;
		}
		else
		{
			sta_end->server[ind] = false;
			sta_end->distance[ind] = MAX_DIS;
		}
	}
	totalcoms = temp_totalcoms;
	min_max(list_v, list_e, sta_end, num_nn, sta_end->total_consume, &temp_cost);
	for (int i = 0; i<num_nn; i++)
	{
		sta_end->solveNode[i] = -1;//solveNode[i]=0;
		if (sta_end->flow[i]>0) sta_end->solveNode[i] = sta_end->flag_black[i];
	}
	cost_simp = simplex(prob, sta_end->solveNode, sta_end->solveLink, sta_end->solveSerRemain);
	check_opt(sta_end, cost_simp);

whole_end:
	return 0;
}
void adj_ansLge(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out)
{
	priority_queue<pri_que> chg_queue;
	pri_que temp_pri;
	Change change;
	unsigned int sum_cst, sum_need, opt;
	short temp_top, oppo_temp, temp_ser_level;
	unsigned short oppo_temp_level, flowOppo, flowTop;

	recv_best(list_v, list_e, sta_end);

	for (ind = 0; ind<num_nn; ind++)
	{
		if (sta_end->flow[ind]>0)
		{
			sum_cst = list_v[ind].base_cost + sta_end->ser_cost[list_v[ind].ser_level];
			sum_need = sta_end->flow[ind];
			temp_pri.id = ind;
			temp_pri.demands = list_v[ind].demand;
			temp_pri.percost = -float(sum_cst) / float(sum_need);
			temp_pri.cost = -float(sum_cst) / float(sum_need);
			chg_queue.push(temp_pri);
		}
	}
	if (jump_time(time_out)) { goto whole_end; }
	while (!chg_queue.empty())
	{
		recv_best(list_v, list_e, sta_end);
		temp_pri = chg_queue.top();
		chg_queue.pop();
		temp_top = temp_pri.id;
		temp_ser_level = list_v[temp_top].ser_level;
		//delete current solution
#ifdef DELETE	
		change.changeNum = 1;
		change.changeIndex[0] = temp_top;
		change.changeOld[0] = temp_ser_level;
		change.changeNew[0] = -1;
		opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
		if (check_opt(sta_end, opt)) {
#ifdef COUTINFO				
			cout << "Delete simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
			cout << "Delete opt" << opt << endl;
#endif
			bkup_best(list_v, list_e, sta_end);
		}
#endif
		change.changeIndex[1] = temp_top;
		change.changeOld[1] = temp_ser_level;
		change.changeNew[1] = -1;
		change.changeNum = 2;
		for (int chg_ind = 0; chg_ind<list_v[temp_top].num_edge; chg_ind++)
		{
			recover_graph(list_v, list_e, sta_end);
			oppo_temp = list_v[temp_top].oppo[chg_ind];
			change.changeIndex[0] = oppo_temp;
			change.changeOld[0] = list_v[oppo_temp].ser_level;
			change.changeNew[1] = -1;
			if (list_v[oppo_temp].ser_level<0 && sta_end->flag_black[oppo_temp] >= 0)
			{
				change.changeNew[0] = temp_ser_level;
#ifdef COUTINFO				
				//cout<<"move start"<<endl; 
#endif
				opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
				if (check_opt(sta_end, opt))
				{
#ifdef COUTINFO				
					cout << "Move simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
					cout << "Move opt" << opt << endl;
#endif
					bkup_best(list_v, list_e, sta_end);
				}
				if (jump_time(time_out)) goto whole_end;
			}
			else if (list_v[oppo_temp].ser_level >= 0 && sta_end->flag_black[oppo_temp] >= 0)
			{
#ifdef MERGE
				oppo_temp_level = list_v[oppo_temp].ser_level;
				flowOppo = sta_end->flow[oppo_temp];
				flowTop = sta_end->flow[temp_top];
				//cout<<"Merge start"<<endl;
				if ((flowOppo + flowTop)>sta_end->ser_upper[sta_end->flag_black[oppo_temp]])
				{
					change.changeNew[0] = sta_end->flag_black[oppo_temp];
					change.changeNew[1] = cal_level(sta_end, (flowOppo + flowTop) - sta_end->ser_upper[sta_end->flag_black[oppo_temp]]);
				}
				else
				{
					change.changeNew[0] = cal_level(sta_end, (flowOppo + flowTop));
					change.changeNew[1] = -1;
				}
				list_v[oppo_temp].ser_level = change.changeNew[0];
				list_v[temp_top].ser_level = change.changeNew[1];
				opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
				if (check_opt(sta_end, opt))
				{
#ifdef COUTINFO				
					cout << "MERGE simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
					cout << "MERGE opt" << opt << endl;
#endif
					bkup_best(list_v, list_e, sta_end);
				}
				if (jump_time(time_out))  goto whole_end;
#endif         
			}
		}
		recv_best(list_v, list_e, sta_end);
#ifdef RELOAD
		short idMaxresi, idMinresi, resiMax, resiMin, resiTemp;
		resiMax = 0; resiMin = sta_end->ser_upper[sta_end->ser_num - 1];
		idMaxresi = -1; idMinresi = -1;
		for (int i = 0; i<num_nn; i++)
		{
			if (sta_end->flow[i])
			{
				resiTemp = sta_end->ser_upper[list_v[i].ser_level] - sta_end->flow[i];
				if (resiTemp<resiMin && (list_v[i].ser_level<sta_end->ser_num - 1))
					//if(resiTemp<resiMin&&list_v[i].ser_level<sta_end->flag_black[i]) 
				{
					idMinresi = i;
					resiMin = resiTemp;
				}
				if (resiTemp>resiMax)
				{
					idMaxresi = i;
					resiMax = resiTemp;
				}
			}
		}
		if (idMaxresi != -1 && idMinresi != -1 && idMinresi != idMaxresi)
		{
			change.changeNum = 2;
			change.changeIndex[0] = idMaxresi;
			change.changeOld[0] = list_v[idMaxresi].ser_level;
			change.changeNew[0] = list_v[idMaxresi].ser_level - 1;
			change.changeIndex[1] = idMinresi;
			change.changeOld[1] = list_v[idMinresi].ser_level;
			change.changeNew[1] = list_v[idMinresi].ser_level + 1;
#ifdef COUTINFO				
			//cout<<"Reload start"<<endl;
#endif
			opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
			if (check_opt(sta_end, opt))
			{
#ifdef COUTINFO				
				cout << "RELOAD simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
				cout << "RELOAD opt" << opt << endl;
#endif
				bkup_best(list_v, list_e, sta_end);
			}
			if (jump_time(time_out)) goto whole_end;
		}
		recv_best(list_v, list_e, sta_end);
		if (jump_time(time_out)) goto whole_end;
#endif
		backup_graph(list_v, list_e, sta_end);
	}
whole_end:;

}

int iter_ansMid2(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out)
{
	unsigned int temp_cost_sum, temp_need_sum;
	short temp_cost, ser_cnt, diff, ind;;
	priority_queue<pri_que> p_quque;
	pri_que temp_pri;
	unsigned long  temp_totalcoms = totalcoms;
	stack<unsigned short> stack_delete;
	unsigned int cost_simp;
	float alpha, alpha_s, alpha_e, alpha_step;
	int temp_ser_level;
	unsigned short output_reduced, index_set = 0, temp_real_ind, stack_top;

	alpha_s = (float)0.5; alpha_e = (float)1.7; alpha_step = (float)0.1;
	alpha = alpha_s;
	for (ind = 0; ind<num_nn; ind++)
	{
		list_v[ind].ser_level = sta_end->flag_black[ind];
		if (sta_end->flag_black[ind] >= 0)
		{
			temp_cost_sum = list_v[ind].base_cost;
			temp_need_sum = list_v[ind].demand;
			temp_ser_level = sta_end->flag_black[ind];
			for (int edge_ind = 0; edge_ind<list_v[ind].num_edge; edge_ind++)
			{
				int edge = list_v[ind].edge[edge_ind];
				temp_cost_sum += (unsigned int)list_e[edge].cost*list_e[edge].limit;
				temp_need_sum += list_e[edge].limit;
			}
			list_v[ind].judge = float(temp_cost_sum + sta_end->ser_cost[temp_ser_level]) / float(temp_need_sum);
		}
		else list_v[ind].judge = MAX_DIS;
	}

	ind = 0;
	while (alpha<alpha_e)
	{
		for (ind = 0; ind<num_link; ind++)
		{
			list_e[ind].cost_dsg[0] = list_e[ind].cost_dsg[1] = list_e[ind].cost;
			list_e[ind].resi_dsg[0] = list_e[ind].resi_dsg[1] = list_e[ind].limit;
			list_e[ind].bnd_dsg[0] = list_e[ind].bnd_dsg[1] = list_e[ind].limit;
		}
		for (ind = 0; ind<num_nn; ind++)
		{
			sta_end->flow[ind] = 0;
			list_v[ind].given = 0;
			if (list_v[ind].ser_level >= 0)
			{
				sta_end->server[ind] = true;
				sta_end->distance[ind] = alpha*list_v[ind].judge;
			}
			else
			{
				sta_end->server[ind] = false;
				sta_end->distance[ind] = MAX_DIS;
			}
		}

		if (jump_time(time_out)) goto whole_end;
		//totalcoms = temp_totalcoms;
		min_max(list_v, list_e, sta_end, num_nn, sta_end->total_consume, &temp_cost);
		for (int i = 0; i<num_nn; i++)
		{
			sta_end->solveNode[i] = -1;//solveNode[i]=0;
			if (sta_end->flow[i]>0) sta_end->solveNode[i] = sta_end->flag_black[i];
		}
		if (jump_time(time_out)) goto whole_end;
		cost_simp = simplex(prob, sta_end->solveNode, sta_end->solveLink, sta_end->solveSerRemain);
		check_opt(sta_end, cost_simp);
		//start the compress
		int output_total = 0;
		for (int i = 0; i<num_nn; i++)
		{
			if (sta_end->solveNode[i] >= 0)
			{
				sta_end->flow[i] = sta_end->ser_upper[sta_end->solveNode[i]] - sta_end->solveSerRemain[i];
				output_total += sta_end->ser_upper[cal_level(sta_end, sta_end->flow[i])];
				temp_pri.id = i;
				temp_pri.percost = float(sta_end->flow[i]) / float(sta_end->ser_upper[cal_level(sta_end, sta_end->flow[i])]);
				temp_pri.demands = cal_level(sta_end, sta_end->flow[i]);
				sta_end->solveNode[i] = temp_pri.demands;
				p_quque.push(temp_pri);
			}
		}

		if (jump_time(time_out)) goto whole_end;
		ind = 0;
		while (!p_quque.empty())
		{
			temp_pri = p_quque.top();
			p_quque.pop();
			sta_end->list_ser[ind] = temp_pri.id;
			ind++;
		}
		ser_cnt = ind;
		output_reduced = output_total;
		index_set = 0;
		for (int ind = 0; ind<ser_cnt; ind++)
		{
			temp_real_ind = sta_end->list_ser[ind];
			index_set++;
			sta_end->solveNode[temp_real_ind]--;
			diff = sta_end->ser_upper[sta_end->solveNode[temp_real_ind]] - sta_end->ser_upper[sta_end->solveNode[temp_real_ind] + 1];
			output_reduced += diff;
			if (output_reduced <= sta_end->total_consume)
				//if(output_reduced<=temp_totalcoms)
			{
				sta_end->solveNode[temp_real_ind]++;
				goto while_chg_end;
			}
			else stack_delete.push(temp_real_ind);
		}
		if (jump_time(time_out)) goto whole_end;

	while_chg_end:
		while (!stack_delete.empty())
		{
			cost_simp = simplex(prob, sta_end->solveNode, sta_end->solveLink, sta_end->solveSerRemain);
			if (cost_simp>1000000)
			{
				stack_top = stack_delete.top();
				stack_delete.pop();
				sta_end->solveNode[stack_top]++;
			}
			else
			{
				while (!stack_delete.empty()) stack_delete.pop();
				goto while_end;
			}
			if (jump_time(time_out)) break;
		}
	while_end:
		check_opt(sta_end, cost_simp);
		alpha = alpha + alpha_step;
		if (jump_time(time_out)) goto whole_end;
	}

whole_end:
	return 0;
}
void adj_ansMid2(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out)
{
	priority_queue<pri_que> chg_queue;
	pri_que temp_pri;
	Change change;
	unsigned int sum_cst, sum_need, opt;
	short temp_top, oppo_temp, temp_ser_level;
	unsigned short oppo_temp_level, flowOppo, flowTop;

	recv_best(list_v, list_e, sta_end);

	for (ind = 0; ind<num_nn; ind++)
	{
		if (sta_end->flow[ind]>0)
		{
			sum_cst = list_v[ind].base_cost + sta_end->ser_cost[list_v[ind].ser_level];
			sum_need = sta_end->flow[ind];
			temp_pri.id = ind;
			temp_pri.demands = list_v[ind].demand;
			temp_pri.percost = -float(sum_cst) / float(sum_need);
			temp_pri.cost = -float(sum_cst) / float(sum_need);
			chg_queue.push(temp_pri);
		}
	}
	if (jump_time(time_out)) { goto whole_end; }
	while (!chg_queue.empty())
	{
		recv_best(list_v, list_e, sta_end);
		temp_pri = chg_queue.top();
		chg_queue.pop();
		temp_top = temp_pri.id;
		temp_ser_level = list_v[temp_top].ser_level;
		//delete current solution
		change.changeNum = 1;
		change.changeIndex[0] = temp_top;
		change.changeOld[0] = temp_ser_level;
		change.changeNew[0] = -1;
		if (jump_time(time_out)) goto whole_end;
		opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
		if (check_opt(sta_end, opt)) {
#ifdef COUTINFO				
			cout << "Delete simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
			cout << "Delete opt" << opt << endl;
#endif
			bkup_best(list_v, list_e, sta_end);
		}
		change.changeIndex[1] = temp_top;
		change.changeOld[1] = temp_ser_level;
		change.changeNew[1] = -1;
		change.changeNum = 2;
		for (int chg_ind = 0; chg_ind<list_v[temp_top].num_edge; chg_ind++)
		{
			recover_graph(list_v, list_e, sta_end);
			oppo_temp = list_v[temp_top].oppo[chg_ind];
			change.changeIndex[0] = oppo_temp;
			change.changeOld[0] = list_v[oppo_temp].ser_level;
			change.changeNew[1] = -1;
			if (list_v[oppo_temp].ser_level<0 && sta_end->flag_black[oppo_temp] >= 0)
			{
				change.changeNew[0] = temp_ser_level;
#ifdef COUTINFO				
				cout << "move start" << endl;
#endif
				opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
				if (check_opt(sta_end, opt))
				{
#ifdef COUTINFO				
					cout << "simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
					cout << "Move opt" << opt << endl;
#endif
					bkup_best(list_v, list_e, sta_end);
				}
				if (jump_time(time_out)) goto whole_end;
			}
			else if (list_v[oppo_temp].ser_level >= 0 && sta_end->flag_black[oppo_temp] >= 0)
			{
#ifdef MERGE
				oppo_temp_level = list_v[oppo_temp].ser_level;
				flowOppo = sta_end->flow[oppo_temp];
				flowTop = sta_end->flow[temp_top];
				//cout<<"Merge start"<<endl;
				if ((flowOppo + flowTop)>sta_end->ser_upper[sta_end->flag_black[oppo_temp]])
				{
					change.changeNew[0] = sta_end->flag_black[oppo_temp];
					change.changeNew[1] = cal_level(sta_end, (flowOppo + flowTop) - sta_end->ser_upper[sta_end->flag_black[oppo_temp]]);
				}
				else
				{
					change.changeNew[0] = cal_level(sta_end, (flowOppo + flowTop));
					change.changeNew[1] = -1;
				}
				list_v[oppo_temp].ser_level = change.changeNew[0];
				list_v[temp_top].ser_level = change.changeNew[1];
				opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
				if (check_opt(sta_end, opt))
				{
#ifdef COUTINFO				
					cout << "simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
					cout << "MERGE opt" << opt << endl;
#endif
					bkup_best(list_v, list_e, sta_end);
				}
				if (jump_time(time_out))  goto whole_end;
#endif         
			}
		}
		recv_best(list_v, list_e, sta_end);
#ifdef RELOAD
		short idMaxresi, idMinresi, resiMax, resiMin, resiTemp;
		resiMax = 0; resiMin = sta_end->ser_upper[sta_end->ser_num - 1];
		idMaxresi = -1; idMinresi = -1;
		for (int i = 0; i<num_nn; i++)
		{
			if (sta_end->flow[i])
			{
				resiTemp = sta_end->ser_upper[list_v[i].ser_level] - sta_end->flow[i];
				if (resiTemp<resiMin && (list_v[i].ser_level<sta_end->ser_num - 1))
					//if(resiTemp<resiMin&&list_v[i].ser_level<sta_end->flag_black[i]) 
				{
					idMinresi = i;
					resiMin = resiTemp;
				}
				if (resiTemp>resiMax)
				{
					idMaxresi = i;
					resiMax = resiTemp;
				}
			}
		}
		if (idMaxresi != -1 && idMinresi != -1 && idMinresi != idMaxresi)
		{
			change.changeNum = 2;
			change.changeIndex[0] = idMaxresi;
			change.changeOld[0] = list_v[idMaxresi].ser_level;
			change.changeNew[0] = list_v[idMaxresi].ser_level - 1;
			change.changeIndex[1] = idMinresi;
			change.changeOld[1] = list_v[idMinresi].ser_level;
			change.changeNew[1] = list_v[idMinresi].ser_level + 1;
#ifdef COUTINFO				
			cout << "Reload start" << endl;
#endif
			opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
			if (check_opt(sta_end, opt))
			{
#ifdef COUTINFO				
				cout << "simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
				cout << "RELOAD opt" << opt << endl;
#endif
				bkup_best(list_v, list_e, sta_end);
			}
			if (jump_time(time_out)) goto whole_end;
		}
		recv_best(list_v, list_e, sta_end);
		if (jump_time(time_out)) goto whole_end;
#endif
		backup_graph(list_v, list_e, sta_end);
	}
whole_end:;

}
int iter_ansLge2(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out)
{
	unsigned int temp_cost_sum, temp_need_sum;
	short temp_cost, ser_cnt, diff, ind;;
	priority_queue<pri_que> p_quque;
	pri_que temp_pri;
	unsigned long  temp_totalcoms = totalcoms;
	stack<unsigned short> stack_delete;
	unsigned int cost_simp;
	float alpha, alpha_s, alpha_e, alpha_step;
	int temp_ser_level;
	unsigned short output_reduced, index_set = 0, temp_real_ind, stack_top;

	alpha_s = (float)0.5; alpha_e = (float)1.7; alpha_step = (float)0.1;
	alpha = alpha_s;
	for (ind = 0; ind<num_nn; ind++)
	{
		list_v[ind].ser_level = sta_end->flag_black[ind];
		if (sta_end->flag_black[ind] >= 0)
		{
			temp_cost_sum = list_v[ind].base_cost;
			temp_need_sum = list_v[ind].demand;
			temp_ser_level = sta_end->flag_black[ind];
			for (int edge_ind = 0; edge_ind<list_v[ind].num_edge; edge_ind++)
			{
				int edge = list_v[ind].edge[edge_ind];
				temp_cost_sum += (unsigned int)list_e[edge].cost*list_e[edge].limit;
				temp_need_sum += list_e[edge].limit;
			}
			list_v[ind].judge = float(temp_cost_sum + sta_end->ser_cost[temp_ser_level]) / float(temp_need_sum);
		}
		else list_v[ind].judge = MAX_DIS;
	}

	ind = 0;

	for (ind = 0; ind<num_link; ind++)
	{
		list_e[ind].cost_dsg[0] = list_e[ind].cost_dsg[1] = list_e[ind].cost;
		list_e[ind].resi_dsg[0] = list_e[ind].resi_dsg[1] = list_e[ind].limit;
		list_e[ind].bnd_dsg[0] = list_e[ind].bnd_dsg[1] = list_e[ind].limit;
	}
	for (ind = 0; ind<num_nn; ind++)
	{
		sta_end->flow[ind] = 0;
		list_v[ind].given = 0;
		if (list_v[ind].ser_level >= 0)
		{
			sta_end->server[ind] = true;
			sta_end->distance[ind] = alpha*list_v[ind].judge;
		}
		else
		{
			sta_end->server[ind] = false;
			sta_end->distance[ind] = MAX_DIS;
		}
	}
	//totalcoms = temp_totalcoms;
	min_max(list_v, list_e, sta_end, num_nn, sta_end->total_consume, &temp_cost);
	for (int i = 0; i<num_nn; i++)
	{
		sta_end->solveNode[i] = -1;//solveNode[i]=0;
		if (sta_end->flow[i]>0) sta_end->solveNode[i] = sta_end->flag_black[i];
	}
	cost_simp = simplex(prob, sta_end->solveNode, sta_end->solveLink, sta_end->solveSerRemain);
	check_opt(sta_end, cost_simp);
	//start the compress
	int output_total = 0;
	for (int i = 0; i<num_nn; i++)
	{
		if (sta_end->solveNode[i] >= 0)
		{
			sta_end->flow[i] = sta_end->ser_upper[sta_end->solveNode[i]] - sta_end->solveSerRemain[i];
			output_total += sta_end->ser_upper[cal_level(sta_end, sta_end->flow[i])];
			temp_pri.id = i;
			temp_pri.percost = float(sta_end->flow[i]) / float(sta_end->ser_upper[cal_level(sta_end, sta_end->flow[i])]);
			temp_pri.demands = cal_level(sta_end, sta_end->flow[i]);
			sta_end->solveNode[i] = temp_pri.demands;
			p_quque.push(temp_pri);
		}
	}
	ind = 0;
	while (!p_quque.empty())
	{
		temp_pri = p_quque.top();
		p_quque.pop();
		sta_end->list_ser[ind] = temp_pri.id;
		ind++;
	}
	ser_cnt = ind;
	output_reduced = output_total;
	index_set = 0;
	for (int ind = 0; ind<ser_cnt; ind++)
	{
		temp_real_ind = sta_end->list_ser[ind];
		index_set++;
		sta_end->solveNode[temp_real_ind]--;
		diff = sta_end->ser_upper[sta_end->solveNode[temp_real_ind]] - sta_end->ser_upper[sta_end->solveNode[temp_real_ind] + 1];
		output_reduced += diff;
		if (output_reduced <= sta_end->total_consume)
			//if(output_reduced<=temp_totalcoms)
		{
			sta_end->solveNode[temp_real_ind]++;
			goto while_chg_end;
		}
		else stack_delete.push(temp_real_ind);
	}
	if (jump_time(time_out)) goto whole_end;

while_chg_end:
	while (!stack_delete.empty())
	{
		cost_simp = simplex(prob, sta_end->solveNode, sta_end->solveLink, sta_end->solveSerRemain);
		if (cost_simp>1000000)
		{
			stack_top = stack_delete.top();
			stack_delete.pop();
			sta_end->solveNode[stack_top]++;
		}
		else
		{
			while (!stack_delete.empty()) stack_delete.pop();
			goto while_end;
		}
		if (jump_time(time_out)) break;
	}
while_end:
	check_opt(sta_end, cost_simp);

whole_end:
	return 0;
}
void adj_ansLge2(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned long totalcoms, unsigned int time_out)
{
	priority_queue<pri_que> chg_queue;
	pri_que temp_pri;
	Change change;
	unsigned int sum_cst, sum_need, opt;
	short temp_top, oppo_temp, temp_ser_level;
	unsigned short oppo_temp_level, flowOppo, flowTop;

	recv_best(list_v, list_e, sta_end);

	for (ind = 0; ind<num_nn; ind++)
	{
		if (sta_end->flow[ind]>0)
		{
			sum_cst = list_v[ind].base_cost + sta_end->ser_cost[list_v[ind].ser_level];
			sum_need = sta_end->flow[ind];
			temp_pri.id = ind;
			temp_pri.demands = list_v[ind].demand;
			temp_pri.percost = -float(sum_cst) / float(sum_need);
			temp_pri.cost = -float(sum_cst) / float(sum_need);
			chg_queue.push(temp_pri);
		}
	}
	if (jump_time(time_out)) { goto whole_end; }
	while (!chg_queue.empty())
	{
		recv_best(list_v, list_e, sta_end);
		temp_pri = chg_queue.top();
		chg_queue.pop();
		temp_top = temp_pri.id;
		temp_ser_level = list_v[temp_top].ser_level;
		//delete current solution
#ifdef DELETE	
		change.changeNum = 1;
		change.changeIndex[0] = temp_top;
		change.changeOld[0] = temp_ser_level;
		change.changeNew[0] = -1;
		opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
		if (check_opt(sta_end, opt)) {
#ifdef COUTINFO				
			cout << "Delete simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
			cout << "Delete opt" << opt << endl;
#endif
			bkup_best(list_v, list_e, sta_end);
		}
#endif
		change.changeIndex[1] = temp_top;
		change.changeOld[1] = temp_ser_level;
		change.changeNew[1] = -1;
		change.changeNum = 2;
		for (int chg_ind = 0; chg_ind<list_v[temp_top].num_edge; chg_ind++)
		{
			recover_graph(list_v, list_e, sta_end);
			oppo_temp = list_v[temp_top].oppo[chg_ind];
			change.changeIndex[0] = oppo_temp;
			change.changeOld[0] = list_v[oppo_temp].ser_level;
			change.changeNew[1] = -1;
			if (list_v[oppo_temp].ser_level<0 && sta_end->flag_black[oppo_temp] >= 0)
			{
				change.changeNew[0] = temp_ser_level;
#ifdef COUTINFO				
				//cout<<"move start"<<endl; 
#endif
				opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
				if (check_opt(sta_end, opt))
				{
#ifdef COUTINFO				
					cout << "Move simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
					cout << "Move opt" << opt << endl;
#endif
					bkup_best(list_v, list_e, sta_end);
				}
				if (jump_time(time_out)) goto whole_end;
			}
			else if (list_v[oppo_temp].ser_level >= 0 && sta_end->flag_black[oppo_temp] >= 0)
			{
#ifdef MERGE
				oppo_temp_level = list_v[oppo_temp].ser_level;
				flowOppo = sta_end->flow[oppo_temp];
				flowTop = sta_end->flow[temp_top];
				//cout<<"Merge start"<<endl;
				if ((flowOppo + flowTop)>sta_end->ser_upper[sta_end->flag_black[oppo_temp]])
				{
					change.changeNew[0] = sta_end->flag_black[oppo_temp];
					change.changeNew[1] = cal_level(sta_end, (flowOppo + flowTop) - sta_end->ser_upper[sta_end->flag_black[oppo_temp]]);
				}
				else
				{
					change.changeNew[0] = cal_level(sta_end, (flowOppo + flowTop));
					change.changeNew[1] = -1;
				}
				list_v[oppo_temp].ser_level = change.changeNew[0];
				list_v[temp_top].ser_level = change.changeNew[1];
				opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
				if (check_opt(sta_end, opt))
				{
#ifdef COUTINFO				
					cout << "MERGE simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
					cout << "MERGE opt" << opt << endl;
#endif
					bkup_best(list_v, list_e, sta_end);
				}
				if (jump_time(time_out))  goto whole_end;
#endif         
			}
		}
		recv_best(list_v, list_e, sta_end);
#ifdef RELOAD
		short idMaxresi, idMinresi, resiMax, resiMin, resiTemp;
		resiMax = 0; resiMin = sta_end->ser_upper[sta_end->ser_num - 1];
		idMaxresi = -1; idMinresi = -1;
		for (int i = 0; i<num_nn; i++)
		{
			if (sta_end->flow[i])
			{
				resiTemp = sta_end->ser_upper[list_v[i].ser_level] - sta_end->flow[i];
				if (resiTemp<resiMin && (list_v[i].ser_level<sta_end->ser_num - 1))
					//if(resiTemp<resiMin&&list_v[i].ser_level<sta_end->flag_black[i]) 
				{
					idMinresi = i;
					resiMin = resiTemp;
				}
				if (resiTemp>resiMax)
				{
					idMaxresi = i;
					resiMax = resiTemp;
				}
			}
		}
		if (idMaxresi != -1 && idMinresi != -1 && idMinresi != idMaxresi)
		{
			change.changeNum = 2;
			change.changeIndex[0] = idMaxresi;
			change.changeOld[0] = list_v[idMaxresi].ser_level;
			change.changeNew[0] = list_v[idMaxresi].ser_level - 1;
			change.changeIndex[1] = idMinresi;
			change.changeOld[1] = list_v[idMinresi].ser_level;
			change.changeNew[1] = list_v[idMinresi].ser_level + 1;
#ifdef COUTINFO				
			//cout<<"Reload start"<<endl;
#endif
			opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
			if (check_opt(sta_end, opt))
			{
#ifdef COUTINFO				
				cout << "RELOAD simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
				cout << "RELOAD opt" << opt << endl;
#endif
				bkup_best(list_v, list_e, sta_end);
			}
			if (jump_time(time_out)) goto whole_end;
		}
		recv_best(list_v, list_e, sta_end);
		if (jump_time(time_out)) goto whole_end;
#endif
		backup_graph(list_v, list_e, sta_end);
	}
whole_end:;

}

int road_milp(pv list_v, pe list_e, pse stt_end, unsigned long totalcoms, char *output, short *solve)
{
	unsigned short cnt = 0;
	stack<unsigned short> reverse;
	string out;
	short temp_edge_f;
	unsigned int time_out = 89000;
	while (totalcoms)
	{
		min_max_q->end = min_max_q->top = 0;
		min_max_q->cntn = 0;
		stt_end->dist = MAX_DIS;
		for (int i = 0; i<num_nn; i++)
		{
			list_v[i].flag_q = false;
			list_v[i].pre = num_nn;
			if (solve[i + 2 * num_link] >= 0 && (stt_end->ser_upper[solve[i + 2 * num_link]] - stt_end->flow[i])>0)
			{
				min_max_q->queue[min_max_q->end] = i;
				min_max_q->end++;
				min_max_q->cntn++;
				list_v[i].flag_q = true;
				//list_v[i].band=MAX_CST;
				if (list_v[i].customer&&list_v[i].given<list_v[i].demand)
				{
					stt_end->dist = 0;
					stt_end->pre = i;
					goto while_end;
				}
			}
		}
		while (nonempty_q(min_max_q))
		{
			temp_pop = min_max_q->queue[min_max_q->top];
			min_max_q->top++;
			min_max_q->cntn--;
			min_max_q->top = min_max_q->top%min_max_q->length;
			//list_v[temp_pop].flag_q=false;

			for (ind = 0; ind<list_v[temp_pop].num_edge; ind++)
			{
				ind_e = list_v[temp_pop].edge[ind];
				temp_edge_f = list_v[temp_pop].edge_flag[ind];
				if (solve[ind_e + temp_edge_f*num_link] == 0) continue;
				temp_visit = list_v[temp_pop].oppo[ind];
				list_v[temp_visit].pre = temp_pop;
				list_v[temp_visit].pre_arc = ind_e;
				list_e[ind_e].path = temp_edge_f;
				if (list_v[temp_visit].customer&&list_v[temp_visit].given<list_v[temp_visit].demand)
				{
					stt_end->dist = 0;
					stt_end->pre = temp_visit;
					goto while_end;
				}
				if (!list_v[temp_visit].flag_q)
				{
					list_v[temp_visit].flag_q = true;
					min_max_q->queue[min_max_q->end] = temp_visit;
					min_max_q->end++;
					min_max_q->cntn++;
					min_max_q->end = min_max_q->end%min_max_q->length;
				}
			}
		}
		endTime = clock(); if (((double)(endTime - startTime) * 1000 / CLOCKS_PER_SEC)>time_out) goto whole_end;
	while_end:
		consume = 0;
		if (stt_end->dist != MAX_DIS)
		{
			temp_ind = stt_end->pre;
			consume = totalcoms;
			if ((list_v[temp_ind].demand - list_v[temp_ind].given)<consume) consume = list_v[temp_ind].demand - list_v[temp_ind].given;
			while (list_v[temp_ind].pre != num_nn)
			{
				temp_arc = list_v[temp_ind].pre_arc;
				temp_p = list_e[temp_arc].path;
				if (consume>solve[temp_arc + num_link*temp_p]) consume = solve[temp_arc + num_link*temp_p];
				temp_ind = list_v[temp_ind].pre;
			}
			//
			unsigned temp_resi = stt_end->ser_upper[solve[temp_ind + 2 * num_link]] - stt_end->flow[temp_ind];
			//cout<<temp_ind<<" server resi "<<temp_resi<<endl;
			if (consume>temp_resi) consume = temp_resi;

			temp_ind = stt_end->pre;
			list_v[temp_ind].given += consume;
			reverse.push(list_v[stt_end->pre].cons_id);
			while (list_v[temp_ind].pre != num_nn)
			{
				reverse.push(temp_ind);
				temp_arc = list_v[temp_ind].pre_arc;
				temp_p = list_e[temp_arc].path;
				list_e[temp_arc].flag_chg = true;
				solve[temp_arc + num_link*temp_p] = solve[temp_arc + num_link*temp_p] - consume;
				if (solve[temp_arc + num_link*temp_p]<0)
				{
					cout << "stop" << endl;
				}
				temp_ind = list_v[temp_ind].pre;
			}
			stt_end->flow[temp_ind] += consume;
			cnt++;
			sprintf(output, "%hd ", temp_ind);
			out = out + output;
			while (!reverse.empty())
			{

				sprintf(output, "%hd ", reverse.top());
				out = out + output;
				reverse.pop();
			}
			totalcoms = totalcoms - consume;
			if (totalcoms) sprintf(output, "%hd %hd \n", consume, stt_end->opt_solve[temp_ind + 2 * num_link]);
			else sprintf(output, "%hd %hd", consume, stt_end->opt_solve[temp_ind + 2 * num_link]);
			out = out + output;
		}

	}
whole_end:
	sprintf(output, "%hd\n\n", cnt);
	string temp_out = output;

	out = temp_out + out;
	strcpy(output, out.c_str());
	/*
	for(int i=0;i<num_nn;i++)
	{
	if(stt_end->opt_solve[i+2*num_link]>=0)
	{
	temp_sum+=stt_end->flow[i];
	cout<<i<<" "<<stt_end->ser_upper[solve[i+2*num_link]]<<" "<<stt_end->flow[i]<<endl;
	}
	if(list_v[i].customer) cout<<list_v[i].demand<<" "<<list_v[i].given<<endl;
	}
	cout<<temp_sum<<endl;*/
	return 0;
}
int min_max(pv list_v, pe list_e, pse stt_end, int total, unsigned long input, short *cost)
{
	short temp_edge_f;
	while (input != 0)
	{
		min_max_q->end = min_max_q->top = 0;
		min_max_q->cntn = 0;
		stt_end->dist = MAX_DIS;
		for (int i = 0; i<total; i++)
		{
			list_v[i].dist = MAX_DIS;
			list_v[i].pre = total;
			//list_v[i].band=0;
			if (stt_end->server[i] && (stt_end->ser_upper[list_v[i].ser_level] - stt_end->flow[i])>0)
			{
				//cout<<i<<" ";
				min_max_q->queue[min_max_q->end] = i;
				min_max_q->end++;
				min_max_q->cntn++;
				min_max_q->top = min_max_q->top%min_max_q->length;
				//cout<<i<<" "<<stt_end->distance[i];
				list_v[i].dist = stt_end->distance[i];

				list_v[i].flag_q = true;
				//list_v[i].band=MAX_CST;
			}
		}
		//cout<<endl;
		while (nonempty_q(min_max_q))
		{
			temp_pop = min_max_q->queue[min_max_q->top];
			min_max_q->top++;
			min_max_q->cntn--;
			min_max_q->top = min_max_q->top%min_max_q->length;
			list_v[temp_pop].flag_q = false;

			if (list_v[temp_pop].customer&&list_v[temp_pop].given != list_v[temp_pop].demand && (stt_end->dist - list_v[temp_pop].dist)>1e-4)
			{
				stt_end->dist = list_v[temp_pop].dist;
				stt_end->pre = temp_pop;
			}

			for (ind = 0; ind<list_v[temp_pop].num_edge; ind++)
			{
				ind_e = list_v[temp_pop].edge[ind];
				temp_edge_f = list_v[temp_pop].edge_flag[ind];
				temp_cost = list_e[ind_e].cost_dsg[temp_edge_f];
				if (temp_cost == MAX_CST) continue;
				//if((list_v[temp_pop].dist+temp_cost)>1.5*stt_end->dist) continue;
				temp_visit = list_v[temp_pop].oppo[ind];

				if (abs((list_v[temp_pop].dist + temp_cost) - list_v[temp_visit].dist)<1e-3) continue;
				if (((list_v[temp_pop].dist + temp_cost)<list_v[temp_visit].dist) && list_v[temp_pop].pre != temp_visit)
				{
					list_v[temp_visit].dist = list_v[temp_pop].dist + temp_cost;
					list_v[temp_visit].pre = temp_pop;
					list_v[temp_visit].pre_arc = ind_e;
					list_e[ind_e].path = temp_edge_f;
					if (!list_v[temp_visit].flag_q)
					{
						list_v[temp_visit].flag_q = true;
						//push_q(min_max_q,temp_visit);
						min_max_q->queue[min_max_q->end] = temp_visit;
						min_max_q->end++;
						min_max_q->cntn++;
						min_max_q->end = min_max_q->end%min_max_q->length;
					}
					if (list_v[temp_visit].customer&&list_v[temp_visit].given != list_v[temp_visit].demand && (stt_end->dist - list_v[temp_visit].dist)>1e-4)
					{
						stt_end->dist = list_v[temp_visit].dist;
						stt_end->pre = temp_visit;
					}
				}
			}
			// while_end: ind=0;
		}

		consume = 0;
		if (stt_end->dist != MAX_DIS)
		{
			temp_ind = stt_end->pre;
			if ((list_v[temp_ind].demand - list_v[temp_ind].given)<input) consume = list_v[temp_ind].demand - list_v[temp_ind].given;
			else consume = input;
			//cout<<endl<<list_v[temp_ind].demand-list_v[temp_ind].given<<"+ ";
			while (list_v[temp_ind].pre != total)
			{
				temp_arc = list_v[temp_ind].pre_arc;
				temp_p = list_e[temp_arc].path;
				if (consume>list_e[temp_arc].bnd_dsg[temp_p]) consume = list_e[temp_arc].bnd_dsg[temp_p];
				temp_ind = list_v[temp_ind].pre;
			}
			unsigned short ser_remain = stt_end->ser_upper[list_v[temp_ind].ser_level] - stt_end->flow[temp_ind];
			//cout<<temp_ind<<" "<<stt_end->flow[temp_ind]<<" "<<list_v[temp_ind].ser_level<<" "<<stt_end->ser_upper[list_v[temp_ind].ser_level]<<endl;
			if (consume>ser_remain) consume = ser_remain;

			temp_ind = stt_end->pre;

			input = input - consume;
			list_v[temp_ind].given += consume;
			while (list_v[temp_ind].pre != total)
			{
				temp_arc = list_v[temp_ind].pre_arc;
				temp_p = list_e[temp_arc].path;
				list_e[temp_arc].flag_chg = true;

				list_e[temp_arc].resi_dsg[temp_p] = list_e[temp_arc].resi_dsg[temp_p] - consume;
				list_e[temp_arc].resi_dsg[1 - temp_p] = list_e[temp_arc].resi_dsg[1 - temp_p] + consume;

				list_e[temp_arc].bnd_dsg[temp_p] = list_e[temp_arc].resi_dsg[temp_p] % list_e[temp_arc].limit;
				list_e[temp_arc].bnd_dsg[1 - temp_p] = list_e[temp_arc].resi_dsg[1 - temp_p] % list_e[temp_arc].limit;

				if (list_e[temp_arc].resi_dsg[temp_p] == list_e[temp_arc].limit)
				{
					list_e[temp_arc].cost_dsg[temp_p] = list_e[temp_arc].cost;
					list_e[temp_arc].cost_dsg[1 - temp_p] = list_e[temp_arc].cost;
					list_e[temp_arc].bnd_dsg[temp_p] = list_e[temp_arc].bnd_dsg[1 - temp_p] = list_e[temp_arc].limit;
				}
				else if (list_e[temp_arc].resi_dsg[temp_p] == 0)
				{
					list_e[temp_arc].cost_dsg[temp_p] = MAX_CST;
					list_e[temp_arc].bnd_dsg[temp_p] = 0;
					list_e[temp_arc].cost_dsg[1 - temp_p] = -list_e[temp_arc].cost;
					list_e[temp_arc].bnd_dsg[1 - temp_p] = list_e[temp_arc].limit;
				}
				else if (list_e[temp_arc].resi_dsg[temp_p]>list_e[temp_arc].limit)
				{
					list_e[temp_arc].cost_dsg[temp_p] = -list_e[temp_arc].cost;
					list_e[temp_arc].cost_dsg[1 - temp_p] = list_e[temp_arc].cost;
				}
				else
				{
					list_e[temp_arc].cost_dsg[temp_p] = list_e[temp_arc].cost;
					list_e[temp_arc].cost_dsg[1 - temp_p] = -list_e[temp_arc].cost;
				}
				temp_ind = list_v[temp_ind].pre;
			}
			stt_end->flow[temp_ind] += consume;
			stt_end->distance[temp_ind] = 0;
			/* if(stt_end->flow[temp_ind] < stt_end->ser_upper[list_v[temp_ind].ser_level])
			{
			stt_end->distance[temp_ind]=0;
			}
			else
			{
			if(list_v[temp_ind].ser_level<(stt_end->ser_num-1))
			{
			stt_end->distance[temp_ind]=0*list_v[temp_ind].judge_distance;
			list_v[temp_ind].ser_level++;
			}
			else
			{
			stt_end->distance[temp_ind]=MAX_DIS;
			}
			}*/
			//cost[0]=consume*stt_end->dist;
		}

	}

	return consume;
}
unsigned int simplex(const Problem *prob, const short *solveNode, unsigned short *solveLink, unsigned short *serverRemaining)
{
#pragma region Initial of Variables
	double min = 0;
	unsigned short numVariable = prob->numLink * 2 + prob->numConNode;
	unsigned short numConstraint = prob->numNetNode;
	double *A = new double[numConstraint * numVariable];
	double *B = new double[numConstraint];
	double *C = new double[numVariable];
	unsigned short *Nindex = new unsigned short[numVariable];
	unsigned short *Bindex = new unsigned short[numConstraint];
	short *R = new short[numVariable + numConstraint];
	for (int i = 0; i < numConstraint; i++)
	{
		if (solveNode[i] >= 0)
		{
			B[i] = prob->nodeNeed[i] + prob->serverAbility[solveNode[i]];
			min += prob->serverCost[solveNode[i]] + prob->nodeCost[i];
		}
		else
			B[i] = prob->nodeNeed[i];
	}
	memset(A, 0, sizeof(double)* numConstraint * numVariable);
	for (int i = 0; i < prob->numLink; i++)
	{
		C[i] = 0 - prob->linkCost[i];
		C[i + prob->numLink] = 0 - prob->linkCost[i];
		if (B[prob->linkIndex1[i]] < 0)
		{
			A[prob->linkIndex1[i] * numVariable + i] = -1;
			A[prob->linkIndex1[i] * numVariable + i + prob->numLink] = 1;
			C[i] -= M;
			C[i + prob->numLink] += M;
		}
		else
		{
			A[prob->linkIndex1[i] * numVariable + i] = 1;
			A[prob->linkIndex1[i] * numVariable + i + prob->numLink] = -1;
			if (B[prob->linkIndex1[i]] >= 0 && B[prob->linkIndex1[i]] != prob->serverAbility[solveNode[prob->linkIndex1[i]]])
			{
				C[i] += M;
				C[i + prob->numLink] -= M;
			}
		}
		if (B[prob->linkIndex2[i]] < 0)
		{
			A[prob->linkIndex2[i] * numVariable + i] = 1;
			A[prob->linkIndex2[i] * numVariable + i + prob->numLink] = -1;
			C[i] += M;
			C[i + prob->numLink] -= M;
		}
		else
		{
			A[prob->linkIndex2[i] * numVariable + i] = -1;
			A[prob->linkIndex2[i] * numVariable + i + prob->numLink] = 1;
			if (B[prob->linkIndex2[i]] >= 0 && B[prob->linkIndex2[i]] != prob->serverAbility[solveNode[prob->linkIndex2[i]]])
			{
				C[i] -= M;
				C[i + prob->numLink] += M;
			}
		}
		Nindex[i] = i;
		Nindex[i + prob->numLink] = i + prob->numLink;
		R[i] = 1;
		R[i + prob->numLink] = 1;
	}
	unsigned short count = 0;
	for (int i = 0; i < numConstraint; i++)
	{
		if (B[i] < 0)
		{
			B[i] = 0 - B[i];
			A[i * numVariable + prob->numLink * 2 + count] = -1;
			C[prob->numLink * 2 + count] = 0 - M;
			Nindex[prob->numLink * 2 + count] = prob->numLink * 2 + count;
			R[prob->numLink * 2 + count] = 1;
			count++;
			min += B[i] * M;
		}
		else if (B[i] > 0 && B[i] != prob->serverAbility[solveNode[i]])
		{
			A[i * numVariable + prob->numLink * 2 + count] = 1;
			C[prob->numLink * 2 + count] = M;
			Nindex[prob->numLink * 2 + count] = prob->numLink * 2 + count;
			R[prob->numLink * 2 + count] = 1;
			count++;
			min += B[i] * M;
		}
		else if (B[i] == 0 && solveNode[i] >= 0)
		{
			A[i * numVariable + prob->numLink * 2 + count] = 1;
			C[prob->numLink * 2 + count] = M;
			Nindex[prob->numLink * 2 + count] = prob->numLink * 2 + count;
			R[prob->numLink * 2 + count] = 1;
			count++;
			min += B[i] * M;
		}
		Bindex[i] = numVariable + i;
		R[numVariable + i] = 0;
	}
#pragma endregion
	unsigned short column;
	unsigned short row;
	double columnNum;
	double rowNum;
	unsigned short tempBN;
	double tempColumn;
	double tempRow;
	double tempc;
	double *tempA = new double[numConstraint];
#pragma region Find Main Column
	columnNum = 1e-10;
	column = M;
	for (int i = 0; i < numVariable; i++)
	{
		if (columnNum < C[i] * R[Nindex[i]])
		{
			columnNum = C[i] * R[Nindex[i]];
			column = i;
		}
	}
#pragma endregion
	while (column != M)
	{
#pragma region Low Bound
		if (R[Nindex[column]] == 1)
		{
			//Find Main Row
			row = M - 1;
			rowNum = M;
			if (Nindex[column] < prob->numLink)
			{
				rowNum = prob->linkBandwidth[Nindex[column]];
				row = M;
			}
			else if (Nindex[column] < prob->numLink * 2)
			{
				rowNum = prob->linkBandwidth[Nindex[column] - prob->numLink];
				row = M;
			}
			else if (Nindex[column] >= prob->numLink * 2 + prob->numConNode &&  Nindex[column] < numVariable)
			{
				rowNum = 1;
				row = M;
			}
			for (int i = 0; i < numConstraint; i++)
			{
				if (A[i * numVariable + column] > 0)
				{
					tempRow = B[i] / A[i * numVariable + column];
					if (tempRow < rowNum)
					{
						rowNum = tempRow;
						row = i;
					}
				}
				else if (A[i * numVariable + column] < 0)
				{
					if (Bindex[i] < prob->numLink)
					{
						tempRow = (B[i] - prob->linkBandwidth[Bindex[i]]) / A[i * numVariable + column];
						if (tempRow < rowNum)
						{
							rowNum = tempRow;
							row = i;
						}
					}
					else if (Bindex[i] < prob->numLink * 2)
					{
						tempRow = (B[i] - prob->linkBandwidth[Bindex[i] - prob->numLink]) / A[i * numVariable + column];
						if (tempRow < rowNum)
						{
							rowNum = tempRow;
							row = i;
						}
					}
					else if (Bindex[i] >= prob->numLink * 2 + prob->numConNode &&  Bindex[i] < numVariable)
					{
						tempRow = (B[i] - 1) / A[i * numVariable + column];
						if (tempRow < rowNum)
						{
							rowNum = tempRow;
							row = i;
						}
					}
				}
			}
			//Don't Need Main Element
			if (row == M)
			{
				for (int i = 0; i < numConstraint; i++)
					B[i] -= A[i * numVariable + column] * rowNum;
				min -= columnNum * rowNum;
				R[Nindex[column]] = -1;
			}
			//Need Main Element
			else
			{
				B[row] -= A[row * numVariable + column] * rowNum;
				if (B[row]<1e-10 &&  B[row]>-1e-10)
					R[Bindex[row]] = 1;
				else
					R[Bindex[row]] = -1;
				R[Nindex[column]] = 0;
				B[row] = rowNum;
				if (A[row * numVariable + column] != 1)
				{
					for (int i = 0; i < column; i++)
						A[row * numVariable + i] /= A[row * numVariable + column];
					for (int i = column + 1; i < numVariable; i++)
						A[row * numVariable + i] /= A[row * numVariable + column];
					tempA[row] = 1 / A[row * numVariable + column];
					A[row * numVariable + column] = 1;
				}
				else
					tempA[row] = 1;
				for (int i = 0; i < row; i++)
				{
					tempColumn = A[i * numVariable + column];
					if (tempColumn != 0)
					{
						for (int j = 0; j < numVariable; j++)
							A[i * numVariable + j] -= tempColumn *  A[row * numVariable + j];
						B[i] -= tempColumn * rowNum;
						tempA[i] = 0 - tempColumn * tempA[row];
					}
					else
						tempA[i] = 0;
				}
				for (int i = row + 1; i < numConstraint; i++)
				{
					tempColumn = A[i * numVariable + column];
					if (tempColumn != 0)
					{
						for (int j = 0; j < numVariable; j++)
							A[i * numVariable + j] -= tempColumn *  A[row * numVariable + j];
						B[i] -= tempColumn * rowNum;
						tempA[i] = 0 - tempColumn * tempA[row];
					}
					else
						tempA[i] = 0;
				}
				tempc = C[column];
				for (int i = 0; i < numVariable; i++)
					C[i] -= tempc *  A[row * numVariable + i];
				min -= columnNum * rowNum;
				for (int i = 0; i < numConstraint; i++)
					A[i * numVariable + column] = tempA[i];
				C[column] = 0 - tempc * tempA[row];
				tempBN = Nindex[column];
				Nindex[column] = Bindex[row];
				Bindex[row] = tempBN;
			}
		}
#pragma endregion
#pragma region Up Bound
		else
		{
			//Find Main Row
			row = M - 1;
			rowNum = M;
			if (Nindex[column] < prob->numLink)
			{
				rowNum = prob->linkBandwidth[Nindex[column]];
				row = M;
			}
			else if (Nindex[column] < prob->numLink * 2)
			{
				rowNum = prob->linkBandwidth[Nindex[column] - prob->numLink];
				row = M;
			}
			else if (Nindex[column] >= prob->numLink * 2 + prob->numConNode &&  Nindex[column] < numVariable)
			{
				rowNum = 1;
				row = M;
			}
			for (int i = 0; i < numConstraint; i++)
			{
				if (A[i * numVariable + column] < 0)
				{
					tempRow = (0 - B[i]) / A[i * numVariable + column];
					if (tempRow < rowNum)
					{
						rowNum = tempRow;
						row = i;
					}
				}
				else if (A[i * numVariable + column] > 0)
				{
					if (Bindex[i] < prob->numLink)
					{
						tempRow = (prob->linkBandwidth[Bindex[i]] - B[i]) / A[i * numVariable + column];
						if (tempRow < rowNum)
						{
							rowNum = tempRow;
							row = i;
						}
					}
					else if (Bindex[i] < prob->numLink * 2)
					{
						tempRow = (prob->linkBandwidth[Bindex[i] - prob->numLink] - B[i]) / A[i * numVariable + column];
						if (tempRow < rowNum)
						{
							rowNum = tempRow;
							row = i;
						}
					}
					else if (Bindex[i] >= prob->numLink * 2 + prob->numConNode &&  Bindex[i] < numVariable)
					{
						tempRow = (1 - B[i]) / A[i * numVariable + column];
						if (tempRow < rowNum)
						{
							rowNum = tempRow;
							row = i;
						}
					}
				}
			}
			//Don't Need Main Element
			if (row == M)
			{
				for (int i = 0; i < numConstraint; i++)
					B[i] += A[i * numVariable + column] * rowNum;
				min -= columnNum * rowNum;
				R[Nindex[column]] = 1;
			}
			//Need Main Element
			else
			{
				B[row] += A[row * numVariable + column] * rowNum;
				if (B[row]<1e-10 &&  B[row]>-1e-10)
					R[Bindex[row]] = 1;
				else
					R[Bindex[row]] = -1;
				R[Nindex[column]] = 0;
				if (A[row * numVariable + column] != 1)
				{
					for (int i = 0; i < column; i++)
						A[row * numVariable + i] /= A[row * numVariable + column];
					for (int i = column + 1; i < numVariable; i++)
						A[row * numVariable + i] /= A[row * numVariable + column];
					tempA[row] = 1 / A[row * numVariable + column];
					A[row * numVariable + column] = 1;
				}
				else
					tempA[row] = 1;
				if (Nindex[column] < prob->numLink)
					B[row] = prob->linkBandwidth[Nindex[column]] - rowNum;
				else if (Nindex[column] < prob->numLink * 2)
					B[row] = prob->linkBandwidth[Nindex[column] - prob->numLink] - rowNum;
				else if (Nindex[column] >= prob->numLink * 2 + prob->numConNode &&  Nindex[column] < numVariable)
					B[row] = 1 - rowNum;
				for (int i = 0; i < row; i++)
				{
					tempColumn = A[i * numVariable + column];
					if (tempColumn != 0)
					{
						for (int j = 0; j < numVariable; j++)
							A[i * numVariable + j] -= tempColumn *  A[row * numVariable + j];
						B[i] += tempColumn * rowNum;
						tempA[i] = 0 - tempColumn * tempA[row];
					}
					else
						tempA[i] = 0;
				}
				for (int i = row + 1; i < numConstraint; i++)
				{
					tempColumn = A[i * numVariable + column];
					if (tempColumn != 0)
					{
						for (int j = 0; j < numVariable; j++)
							A[i * numVariable + j] -= tempColumn *  A[row * numVariable + j];
						B[i] += tempColumn * rowNum;
						tempA[i] = 0 - tempColumn * tempA[row];
					}
					else
						tempA[i] = 0;
				}
				tempc = C[column];
				for (int i = 0; i < numVariable; i++)
					C[i] -= tempc *  A[row * numVariable + i];
				min -= columnNum * rowNum;
				for (int i = 0; i < numConstraint; i++)
					A[i * numVariable + column] = tempA[i];
				C[column] = 0 - tempc * tempA[row];
				tempBN = Nindex[column];
				Nindex[column] = Bindex[row];
				Bindex[row] = tempBN;
			}
		}
#pragma endregion
#pragma region Find Main Column
		columnNum = 0;
		column = M;
		for (int i = 0; i < numVariable; i++)
		{
			if (columnNum < C[i] * R[Nindex[i]])
			{
				columnNum = C[i] * R[Nindex[i]];
				column = i;
			}
		}
#pragma endregion
	}
#pragma region Variable Assignment
	memset(solveLink, 0, sizeof(short)* prob->numLink * 2);
	memset(serverRemaining, 0, sizeof(short)* prob->numNetNode);
	for (int i = 0; i < numConstraint; i++)
	{
		if (B[i] != 0)
		{
			if (Bindex[i] < prob->numLink * 2)
				solveLink[Bindex[i]] = (short)(B[i] + 0.5);
			else
				serverRemaining[i] = B[i];
			if (Bindex[i] >= numVariable)
			{
				if (prob->nodeNeed[Bindex[i] - numVariable] != 0)
					min = 2000000000;
			}
		}
	}
	for (int i = 0; i < numVariable; i++)
	{
		if (Nindex[i] < prob->numLink)
			solveLink[Nindex[i]] = R[Nindex[i]] == 1 ? 0 : prob->linkBandwidth[Nindex[i]];
		else if (Nindex[i] < prob->numLink * 2)
			solveLink[Nindex[i]] = R[Nindex[i]] == 1 ? 0 : prob->linkBandwidth[Nindex[i] - prob->numLink];
	}
#pragma endregion
	delete[] A;
	delete[] B;
	delete[] C;
	delete[] Nindex;
	delete[] Bindex;
	delete[] R;
	delete[] tempA;
	return (int)(min + 0.5);
}
unsigned int simplex_simple(const Problem *prob, const short *solveNode)
{
	unsigned short *solveLink = (unsigned short*)malloc(2 * num_link*sizeof(unsigned short));
	unsigned short *serverRemaining = (unsigned short*)malloc(num_nn*sizeof(unsigned short));
	return simplex(prob, solveNode, solveLink, serverRemaining);
}
void DelRelAdj(pv list_v, pe list_e, pse sta_end, Problem *prob, unsigned int time_out)
{

	priority_queue<pri_que> chg_queue;
	queue<unsigned short> QSer;
	queue<unsigned short> QNeiFir;
	unsigned short NeiSer[100];
	pri_que temp_pri;
	Change change;
	unsigned int sum_cst, sum_need, opt;
	short temp_top, oppo_temp, temp_ser_level, flowDel, flowDelOld, difDel, indIteDel;
	unsigned short oppo_temp_level, flowOppo, flowTop, indexSer, numSearch;

	recv_best(list_v, list_e, sta_end);
	backup_graph(list_v, list_e, sta_end);
	for (ind = 0; ind<num_nn; ind++){
		if (sta_end->flow[ind]>0){
			sum_cst = list_v[ind].base_cost + sta_end->ser_cost[list_v[ind].ser_level];
			sum_need = sta_end->flow[ind];
			temp_pri.id = ind;
			temp_pri.demands = list_v[ind].demand;
			temp_pri.percost = -float(sum_cst) / float(sum_need);
			temp_pri.cost = -float(sum_cst) / float(sum_need);
			chg_queue.push(temp_pri);
		}
	}
	if (jump_time(time_out)) { goto whole_end; }

	while (!chg_queue.empty())
	{
		recv_best(list_v, list_e, sta_end);
		backup_graph(list_v, list_e, sta_end);
		temp_pri = chg_queue.top();
		chg_queue.pop();
		temp_top = temp_pri.id;
		temp_ser_level = list_v[temp_top].ser_level;

		//delete current solution
		change.changeNum = 2;
		change.changeIndex[1] = temp_top;
		change.changeOld[1] = temp_ser_level;
		change.changeNew[1] = -1;

		list_v[temp_top].flag_q = true;

		for (int chg_ind = 0; chg_ind<list_v[temp_top].num_edge; chg_ind++){
			oppo_temp = list_v[temp_top].oppo[chg_ind];
			QNeiFir.push(oppo_temp);
			list_v[oppo_temp].flag_q = true;
		}
		indexSer = 0;
		while (!QNeiFir.empty())
		{
			short tempSec = QNeiFir.front();
			QNeiFir.pop();
			for (int chg_ind = 0; chg_ind<list_v[tempSec].num_edge; chg_ind++)
			{
				oppo_temp = list_v[tempSec].oppo[chg_ind];
				if (list_v[oppo_temp].ser_level >= 0 && !list_v[oppo_temp].flag_q)
				{
					//cout<<oppo_temp<<" ";
					NeiSer[indexSer] = oppo_temp;
					indexSer++;
					list_v[oppo_temp].flag_q = true;
				}
			}
		}
		if (jump_time(time_out)) goto whole_end;
		for (int i = 0; i<num_nn; i++) list_v[i].flag_q = false;
		//cout << "indexSer is " << indexSer << endl;
		for (int indser = 0; indser<indexSer; indser++)
		{
			recover_graph(list_v, list_e, sta_end);
			oppo_temp = NeiSer[indser];
			change.changeIndex[0] = oppo_temp;
			change.changeOld[0] = list_v[oppo_temp].ser_level;
			change.changeNew[1] = -1;
			if (list_v[oppo_temp].ser_level<0 && sta_end->flag_black[oppo_temp] >= 0)
			{
				change.changeNew[0] = temp_ser_level;
#ifdef COUTINFO				
				cout << "move start" << endl;
#endif
				opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
				if (check_opt(sta_end, opt))
				{
#ifdef COUTINFO				
					cout << "simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
					cout << "Move opt" << opt << endl;
#endif
					//cout << "simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
					//cout << "Move opt" << opt << endl;
					bkup_best(list_v, list_e, sta_end);
				}
				if (jump_time(time_out)) goto whole_end;
			}
			else if (list_v[oppo_temp].ser_level >= 0 && sta_end->flag_black[oppo_temp] >= 0)
			{
#ifdef MERGE
				oppo_temp_level = list_v[oppo_temp].ser_level;
				flowOppo = sta_end->flow[oppo_temp];
				flowTop = sta_end->flow[temp_top];
				//cout<<"Merge start"<<endl;
				if ((flowOppo + flowTop)>sta_end->ser_upper[sta_end->flag_black[oppo_temp]])
				{
					change.changeNew[0] = sta_end->flag_black[oppo_temp];
					change.changeNew[1] = cal_level(sta_end, (flowOppo + flowTop) - sta_end->ser_upper[sta_end->flag_black[oppo_temp]]);
				}
				else
				{
					change.changeNew[0] = cal_level(sta_end, (flowOppo + flowTop));
					change.changeNew[1] = -1;
				}
				list_v[oppo_temp].ser_level = change.changeNew[0];
				list_v[temp_top].ser_level = change.changeNew[1];
				opt = resi_map(list_v, list_e, sta_end, sta_end->solveLink, &change);
				if (check_opt(sta_end, opt))
				{
#ifdef COUTINFO				
					cout << "simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
					cout << "MERGE opt" << opt << endl;
#endif
					//cout << "simplex" << simplex_simple(prob, sta_end->solveNode) << endl;
					//cout << "MERGE opt" << opt << endl;
					bkup_best(list_v, list_e, sta_end);
				}
				if (jump_time(time_out))  goto whole_end;
#endif         
			}
		}
	while_end:;
	}
whole_end:;

}
void cnvtGraph(pv list_v, pe list_e, pse sta_end)
{
	for (int ind_ver = 0; ind_ver<num_nn; ind_ver++)
	{
		list_v[ind_ver].given = list_v[ind_ver].demand;
		list_v[ind_ver].ser_level = sta_end->opt_solve[ind_ver + 2 * num_link];
		sta_end->flow[ind_ver] = 0;
		if (list_v[ind_ver].ser_level >= 0)
		{
			sta_end->flow[ind_ver] = list_v[ind_ver].demand;
		}
	}
	unsigned short temp_s, temp_e;
	for (ind = 0; ind<num_link; ind++)
	{
		if (sta_end->opt_solve[ind]>0)
		{
			list_e[ind].resi_dsg[0] = list_e[ind].limit - sta_end->opt_solve[ind];
			list_e[ind].resi_dsg[1] = list_e[ind].limit + sta_end->opt_solve[ind];
			list_e[ind].bnd_dsg[0] = list_e[ind].resi_dsg[0] % list_e[ind].limit;
			list_e[ind].bnd_dsg[1] = list_e[ind].resi_dsg[1] % list_e[ind].limit;

			temp_s = list_e[ind].ind_s;
			temp_e = list_e[ind].ind_e;
			if (list_v[temp_s].ser_level >= 0) sta_end->flow[temp_s] += sta_end->opt_solve[ind];
			if (list_v[temp_e].ser_level >= 0) sta_end->flow[temp_e] -= sta_end->opt_solve[ind];

			if (list_e[ind].resi_dsg[0]<list_e[ind].limit&&list_e[ind].resi_dsg[0]>0)
			{
				list_e[ind].cost_dsg[0] = list_e[ind].cost;
				list_e[ind].cost_dsg[1] = -list_e[ind].cost;
			}
			else if (list_e[ind].resi_dsg[0] == 0)
			{
				list_e[ind].cost_dsg[0] = MAX_CST;
				list_e[ind].bnd_dsg[0] = 0;
				list_e[ind].cost_dsg[1] = -list_e[ind].cost;
				list_e[ind].bnd_dsg[1] = list_e[ind].limit;
			}
		}
		else if (sta_end->opt_solve[ind + num_link]>0)
		{
			list_e[ind].resi_dsg[1] = list_e[ind].limit - sta_end->opt_solve[ind + num_link];
			list_e[ind].resi_dsg[0] = list_e[ind].limit + sta_end->opt_solve[ind + num_link];
			list_e[ind].bnd_dsg[0] = list_e[ind].resi_dsg[0] % list_e[ind].limit;
			list_e[ind].bnd_dsg[1] = list_e[ind].resi_dsg[1] % list_e[ind].limit;

			temp_s = list_e[ind].ind_s;
			temp_e = list_e[ind].ind_e;
			if (list_v[temp_s].ser_level >= 0) sta_end->flow[temp_s] -= sta_end->opt_solve[ind + num_link];
			if (list_v[temp_e].ser_level >= 0) sta_end->flow[temp_e] += sta_end->opt_solve[ind + num_link];

			if (list_e[ind].resi_dsg[1]<list_e[ind].limit&&list_e[ind].resi_dsg[1]>0)
			{
				list_e[ind].cost_dsg[1] = list_e[ind].cost;
				list_e[ind].cost_dsg[0] = -list_e[ind].cost;
			}
			else if (list_e[ind].resi_dsg[1] == 0)
			{
				list_e[ind].cost_dsg[1] = MAX_CST;
				list_e[ind].bnd_dsg[1] = 0;
				list_e[ind].cost_dsg[0] = -list_e[ind].cost;
				list_e[ind].bnd_dsg[0] = list_e[ind].limit;
			}
		}
		else if (sta_end->opt_solve[ind + num_link] == 0 && sta_end->opt_solve[ind] == 0)
		{
			list_e[ind].cost_dsg[0] = list_e[ind].cost_dsg[1] = list_e[ind].cost;
			list_e[ind].resi_dsg[0] = list_e[ind].resi_dsg[1] = list_e[ind].limit;
			list_e[ind].bnd_dsg[0] = list_e[ind].bnd_dsg[1] = list_e[ind].limit;
		}
	}
	backup_graph(list_v, list_e, sta_end);
	bkup_best(list_v, list_e, sta_end);


}
unsigned int resi_map(pv list_v, pe list_e, pse sta_end, unsigned short *solveLink, const Change *change)
{
	clock_t startTemp=clock();
	cntResi++;
	short index_change, level_old, level_new, flow_diff;
	unsigned int input = 0;
	short temp_edge_f;
	int cost_sum;
	unsigned short ser_remain;

	for (int i = 0; i<change->changeNum; i++)
	{
		index_change = change->changeIndex[i];
		level_old = change->changeOld[i];
		level_new = change->changeNew[i];
		//cout<<index_change<<" "<<level_old<<"->"<<level_new<<endl;
		//cal the flow diff
		if (level_new >= level_old) flow_diff = 0;
		else if (level_new >= 0) flow_diff = sta_end->flow[index_change] - sta_end->ser_upper[level_new];
		else flow_diff = sta_end->flow[index_change];
		//change the graph
		list_v[index_change].given -= flow_diff;
		sta_end->flow[index_change] -= flow_diff;
		list_v[index_change].ser_level = change->changeNew[i];
		sta_end->server[index_change] = 1;
		//cal the input
		input += flow_diff;

	}
	//cout<<"input is "<<input<<endl;
	while (input != 0)
	{
		min_max_q->end = min_max_q->top = 0;
		min_max_q->cntn = 0;
		sta_end->dist = MAX_DIS;
		for (int i = 0; i<num_nn; i++)
		{
			list_v[i].dist = MAX_DIS;
			list_v[i].pre = num_nn;
			if (list_v[i].ser_level == -1 || (sta_end->ser_upper[list_v[i].ser_level] - sta_end->flow[i]) == 0) continue;
			min_max_q->queue[min_max_q->end] = i;
			min_max_q->end++;
			min_max_q->cntn++;
			min_max_q->top = min_max_q->top%min_max_q->length;

			list_v[i].dist = 0;
			list_v[i].flag_q = true;
		}
		while (nonempty_q(min_max_q))
		{
			temp_pop = min_max_q->queue[min_max_q->top];
			min_max_q->top++;
			min_max_q->cntn--;
			min_max_q->top = min_max_q->top%min_max_q->length;
			list_v[temp_pop].flag_q = false;

			if (list_v[temp_pop].given != list_v[temp_pop].demand && (sta_end->dist - list_v[temp_pop].dist)>0)
			{
				sta_end->dist = list_v[temp_pop].dist;
				sta_end->pre = temp_pop;
			}

			for (ind = 0; ind<list_v[temp_pop].num_edge; ind++)
			{
				ind_e = list_v[temp_pop].edge[ind];
				temp_edge_f = list_v[temp_pop].edge_flag[ind];
				if ((int)list_e[ind_e].cost_dsg[temp_edge_f] == MAX_CST) continue;
				temp_visit = list_v[temp_pop].oppo[ind];
				temp_cost = list_e[ind_e].cost_dsg[temp_edge_f];
				if ((((int)list_v[temp_pop].dist + (int)temp_cost)<(int)list_v[temp_visit].dist))
				{
					list_v[temp_visit].dist = list_v[temp_pop].dist + temp_cost;
					list_v[temp_visit].pre = temp_pop;
					list_v[temp_visit].pre_arc = ind_e;
					list_e[ind_e].path = temp_edge_f;
					if (!list_v[temp_visit].flag_q)
					{
						list_v[temp_visit].flag_q = true;
						min_max_q->queue[min_max_q->end] = temp_visit;
						min_max_q->end++;
						min_max_q->cntn++;
						min_max_q->end = min_max_q->end%min_max_q->length;
					}
				}
			}
		}
		consume = 0;
		if (sta_end->dist != MAX_DIS)
		{
			temp_ind = sta_end->pre;
			consume = input;
			if ((list_v[temp_ind].demand - list_v[temp_ind].given)<input) consume = (unsigned int)(list_v[temp_ind].demand - list_v[temp_ind].given);
			while (list_v[temp_ind].pre != num_nn)
			{
				temp_arc = list_v[temp_ind].pre_arc;
				temp_p = list_e[temp_arc].path;
				if (consume>list_e[temp_arc].bnd_dsg[temp_p]) consume = list_e[temp_arc].bnd_dsg[temp_p];
				temp_ind = list_v[temp_ind].pre;
			}
			ser_remain = sta_end->ser_upper[list_v[temp_ind].ser_level] - sta_end->flow[temp_ind];
			if (consume>ser_remain) consume = ser_remain;
			temp_ind = sta_end->pre;
			input = input - consume;
			list_v[temp_ind].given += consume;
			while (list_v[temp_ind].pre != num_nn)
			{
				temp_arc = list_v[temp_ind].pre_arc;
				temp_p = list_e[temp_arc].path;

				list_e[temp_arc].resi_dsg[temp_p] = list_e[temp_arc].resi_dsg[temp_p] - consume;
				list_e[temp_arc].resi_dsg[1 - temp_p] = list_e[temp_arc].resi_dsg[1 - temp_p] + consume;

				list_e[temp_arc].bnd_dsg[temp_p] = list_e[temp_arc].resi_dsg[temp_p] % list_e[temp_arc].limit;
				list_e[temp_arc].bnd_dsg[1 - temp_p] = list_e[temp_arc].resi_dsg[1 - temp_p] % list_e[temp_arc].limit;

				if (list_e[temp_arc].resi_dsg[temp_p]>list_e[temp_arc].limit)
				{
					list_e[temp_arc].cost_dsg[temp_p] = -list_e[temp_arc].cost;
					list_e[temp_arc].cost_dsg[1 - temp_p] = list_e[temp_arc].cost;
				}
				else if (list_e[temp_arc].resi_dsg[temp_p]<list_e[temp_arc].limit&&list_e[temp_arc].resi_dsg[temp_p]>0)
				{
					list_e[temp_arc].cost_dsg[temp_p] = list_e[temp_arc].cost;
					list_e[temp_arc].cost_dsg[1 - temp_p] = -list_e[temp_arc].cost;
				}
				else if (list_e[temp_arc].resi_dsg[temp_p] == list_e[temp_arc].limit)
				{
					list_e[temp_arc].cost_dsg[temp_p] = list_e[temp_arc].cost_dsg[1 - temp_p] = list_e[temp_arc].cost;
					list_e[temp_arc].bnd_dsg[temp_p] = list_e[temp_arc].bnd_dsg[1 - temp_p] = list_e[temp_arc].limit;
				}
				else
				{
					list_e[temp_arc].cost_dsg[temp_p] = MAX_CST;
					//list_e[temp_arc].bnd_dsg[temp_p]=0;
					list_e[temp_arc].cost_dsg[1 - temp_p] = -list_e[temp_arc].cost;
					list_e[temp_arc].bnd_dsg[1 - temp_p] = list_e[temp_arc].limit;
				}
				temp_ind = list_v[temp_ind].pre;
			}
			//coutNodeflow( sta_end, list_v);
			sta_end->flow[temp_ind] += consume;
			if (sta_end->flow[temp_ind] < sta_end->ser_upper[list_v[temp_ind].ser_level]) sta_end->distance[temp_ind] = 0;
			else sta_end->distance[temp_ind] = MAX_DIS;
		}
		if (consume == 0)
		{
			return 2000000000;
		}
	}
	cost_sum = 0;
	for (int ind = 0; ind<num_link; ind++)
	{
		cost_sum += abs(list_e[ind].resi_dsg[1] - list_e[ind].resi_dsg[0]) / 2 * list_e[ind].cost;
		if ((list_e[ind].limit - list_e[ind].resi_dsg[0])>0)
		{
			sta_end->solveLink[ind] = list_e[ind].limit - list_e[ind].resi_dsg[0];
			sta_end->solveLink[ind + num_link] = 0;
		}
		else
		{
			sta_end->solveLink[ind] = 0;
			sta_end->solveLink[ind + num_link] = list_e[ind].resi_dsg[0] - list_e[ind].limit;
		}
	}
	for (int i = 0; i<num_nn; i++)
	{
		sta_end->solveNode[i] = -1;
		if (sta_end->flow[i]>0)
		{
			cost_sum += list_v[i].base_cost + sta_end->ser_cost[list_v[i].ser_level];
			//sta_end->solveNode[i]=list_v[i].ser_level;
			sta_end->solveNode[i] = cal_level(sta_end, sta_end->flow[i]);
		}
	}
	clock_t endTemp=clock();
	cout<<(double)(endTemp - startTemp) * 1000 / CLOCKS_PER_SEC<<"ms"<<endl;
	return cost_sum;
}
inline void backup_graph(pv list_v, pe list_e, pse sta_end)
{
	//backup the current graph before the resi_map
	memcpy(sta_end->flow_bkup, sta_end->flow, num_nn*sizeof(unsigned int));
	for (int i = 0; i<num_link; i++)
	{
		list_e[i].cost_dsg[2] = list_e[i].cost_dsg[0]; list_e[i].cost_dsg[3] = list_e[i].cost_dsg[1];
		list_e[i].resi_dsg[2] = list_e[i].resi_dsg[0]; list_e[i].resi_dsg[3] = list_e[i].resi_dsg[1];
		list_e[i].bnd_dsg[2] = list_e[i].bnd_dsg[0]; list_e[i].bnd_dsg[3] = list_e[i].bnd_dsg[1];
	}
	for (int i = 0; i<num_nn; i++)
	{
		sta_end->serLevbkup[i] = list_v[i].ser_level;
		sta_end->given_bkup[i] = list_v[i].given;
	}
}
inline void recover_graph(pv list_v, pe list_e, pse sta_end)
{
	memcpy(sta_end->flow, sta_end->flow_bkup, num_nn*sizeof(unsigned int));
	//memccpy(sta_end->flow,sta_end->flow_bkup,0,num_nn*sizeof(unsigned int));
	for (int i = 0; i<num_link; i++)
	{
		list_e[i].cost_dsg[0] = list_e[i].cost_dsg[2]; list_e[i].cost_dsg[1] = list_e[i].cost_dsg[3];
		list_e[i].resi_dsg[0] = list_e[i].resi_dsg[2]; list_e[i].resi_dsg[1] = list_e[i].resi_dsg[3];
		list_e[i].bnd_dsg[0] = list_e[i].bnd_dsg[2]; list_e[i].bnd_dsg[1] = list_e[i].bnd_dsg[3];
	}
	for (int i = 0; i<num_nn; i++)
	{
		list_v[i].given = sta_end->given_bkup[i];
		list_v[i].ser_level = sta_end->serLevbkup[i];
	}
}
inline void bkup_best(pv list_v, pe list_e, pse sta_end)
{
	//backup the current graph before the resi_map
	memcpy(sta_end->flow_best, sta_end->flow, num_nn*sizeof(unsigned int));
	for (int i = 0; i<num_link; i++)
	{
		list_e[i].cost_dsg[4] = list_e[i].cost_dsg[0]; list_e[i].cost_dsg[5] = list_e[i].cost_dsg[1];
		list_e[i].resi_dsg[4] = list_e[i].resi_dsg[0]; list_e[i].resi_dsg[5] = list_e[i].resi_dsg[1];
		list_e[i].bnd_dsg[4] = list_e[i].bnd_dsg[0]; list_e[i].bnd_dsg[5] = list_e[i].bnd_dsg[1];
	}
	for (int i = 0; i<num_nn; i++)
	{
		sta_end->given_best[i] = list_v[i].given;
		sta_end->serLevbest[i] = cal_level(sta_end, sta_end->flow[i]);
	}
}
inline void recv_best(pv list_v, pe list_e, pse sta_end)
{
	memcpy(sta_end->flow, sta_end->flow_best, num_nn*sizeof(unsigned int));
	for (int i = 0; i<num_link; i++)
	{
		list_e[i].cost_dsg[0] = list_e[i].cost_dsg[4]; list_e[i].cost_dsg[1] = list_e[i].cost_dsg[5];
		list_e[i].resi_dsg[0] = list_e[i].resi_dsg[4]; list_e[i].resi_dsg[1] = list_e[i].resi_dsg[5];
		list_e[i].bnd_dsg[0] = list_e[i].bnd_dsg[4]; list_e[i].bnd_dsg[1] = list_e[i].bnd_dsg[5];
	}
	for (int i = 0; i<num_nn; i++)
	{
		list_v[i].given = sta_end->given_best[i];
		list_v[i].ser_level = sta_end->serLevbest[i];
	}
}
int cal_cost(pse sta_end, unsigned short ser_flow)
{
	if (ser_flow == 0) return 0;
	for (int i = 0; i<sta_end->ser_num; i++)
	{
		if (ser_flow <= sta_end->ser_upper[i]) return sta_end->ser_cost[i];
	}
	return sta_end->ser_cost[sta_end->ser_num - 1];
}
int cal_level(pse sta_end, unsigned short ser_flow)
{
	if (ser_flow == 0) return -1;
	for (int i = 0; i<sta_end->ser_num; i++)
	{
		if (ser_flow <= sta_end->ser_upper[i]) return i;
	}
	return sta_end->ser_num - 1;
}


