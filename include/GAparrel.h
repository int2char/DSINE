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

class NewGAParrel
{
public:
	taskPath* PathSets;
	int taskd;
	Graph &G;
	vector<service>serv;
	int T,E,W,M;
	double *x,*y,*u,*f;
	double *dev_x,*dev_y,*dev_u,*dev_f;
	double*sum,*dev_sum;
	int *paths,*dev_paths;
	int F;
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
		cout<<"t is "<<T<<endl;
		sum=new double[T/512+1];
		E=G.m;
		W=4;
		PathSets=new taskPath[T];
		int p=0;
		for(int i=0;i<ser.size();i++)
			{
				if(flag[i]>0)
					{
						PathSets[p++]=_PathSets[i];
						for(int j=0;j<W;j++)
						{
							int k=0;
							while(PathSets[p-1].Pathset[j][k]>0)k++;
							if(k>M)M=k;
							if(k==0)cout<<"asdasdasd "<<i<<endl;
						}
				}
			}

		paths=new int[M*W*T];
		cout<<"over"<<endl;
		for(int i=0;i<T;i++)
		{
			int ioff=i*W*M;
			for(int j=0;j<W;j++)
			{
				int joff=j*M;
				int k=0;
				while(PathSets[i].Pathset[j][k]>0){
					paths[ioff+joff+k]=PathSets[i].Pathset[j][k];
					k++;
				}
				while(k<M)
				{
					paths[ioff+joff+k]=-1;
					k++;
				}
				if(PathSets[i].Pathset[j][0]<0)cout<<"erro!!!!!!"<<endl;
			}
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
		cout<<"hahaha"<<endl;
	};
public:
	void Cudamalloc();
	vector<pair<string,float> > GAsearch();
};
#endif