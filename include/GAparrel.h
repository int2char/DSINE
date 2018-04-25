#ifndef GAPARREL
#define GAPARREL
#include"Graph.h"
#include "service.h"
#include"taskPath.h"
#include"valuemark.h"
#include"algorithm"
#include"BFS.h"
#include"PathArrange.h"
#include"time.h"
#include"float.h"
struct Ldge{
	int head,tail;
	Ldge(){};
	Ldge(int h,int t):head(h),tail(t){};
};
class NewGAParrel
{
public:
	std::vector<vector<pair<int,int>>>gugu;
	std::vector<set<int>> sset;
	Graph &G;
	vector<service>serv;
	int T,E,W,M,NN;
	double *x,*y,*u,*f,*dev_f;
	double *dev_x,*dev_y,*dev_u,*d,*dev_d;
	double*sum,*dev_sum;
	int *paths,*dev_paths;
	int F;
	int *p,*dev_p;
	Ldge *edge,*dev_edge;
	int *dev_m,*m;
public:
	NewGAParrel(vector<service>&ser, taskPath*_PathSets, Graph &_G):G(_G)
	{
		std::vector<service> nser;
		int c=0;
		std::vector<int> flag(ser.size(),-1);
		for(int i=0;i<ser.size();i++)
			{	
				int co=0;
				if(_PathSets[i].num>=4)
					while(_PathSets[i].Pathset[co][0]>=0&&co<4)co++;
				if(i!=77&&i!=168&&co>=4)
					{	
						flag[i]=1;
						serv.push_back(ser[i]),c++;
					}
			}
		T=c;
		sum=new double[T/512+1];
		E=G.m;
		W=4;
		NN=G.n;
		std::vector<set<int>> tmps(G.n,set<int>());
		sset=tmps;
		std::vector<vector<pair<int,int>>> gu(G.n,vector<pair<int,int>>());
		gugu=gu;
		for(int i=0;i<ser.size();i++)
		{
			sset[ser[i].s].insert(ser[i].t);
			gugu[ser[i].s].push_back(make_pair(ser[i].t,i));

		}
		x=new double[T+1];
		y=new double[T];
		for(int i=0;i<T;i++)
		{
			y[i]=100;
			x[i]=100;
		}
		x[T]=DBL_MAX;
		u=new double[E];
		for(int i=0;i<E;i++)
			u[i]=0.1;
		f=new double[E];
		d=new double[NN*NN];
		p=new int[NN*NN];
		m=new int;
		*m=0;
		Ldge* edge=new Ldge[E];
		for(int i=0;i<E;i++)
		{
			edge[i].head=G.incL[i].head;
			edge[i].tail=G.incL[i].tail;
		}
		
	};
public:
	void Cudamalloc();
	vector<pair<string,float> > GAsearch();
};
#endif