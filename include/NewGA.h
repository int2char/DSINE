/*
 * this is a header and implement of serial GA
 * no correlate .cpp file exit;
 * */
#include"Graph.h"
#include "service.h"
#include"taskPath.h"
#include"valuemark.h"
#include"curand_kernel.h"
#include"iostream"
#include <fstream>
#include"const.h"
#include<algorithm>
#include<vector>
#include<math.h>
#include"BFS.h"
#include"time.h"
#include"PathArrange.h"
#include"float.h"
#include<map>
using namespace std;
#define power .5
#define NNN 10000
double distan[NNN];
int pedge[NNN];
struct FitVM {
	double value;
	int mark;
	FitVM(double _value=0,double _mark=0){
		value = _value;
		mark = _mark;
	}
};
class NewGA
{
	private:
	taskPath* PathSets;
	Graph &G;
	int*st;
	int*te;
	double*demand;
	int**chormes;
	int**childs;
	int**monsters;
	double*capacity;
	double best;
	int totalweight;
	vector<FitVM> Fits;
	vector<pair<int,vector<int>>> Result;
	vector<pair<double,int>>answers;
	public:
		NewGA( Graph &_G):G(_G){
			st = (int*)malloc(Task*sizeof(int));
			te = (int*)malloc(Task*sizeof(int));
			demand = (double*)malloc(Task*sizeof(double));
			capacity = (double*)malloc(sizeof(double)*EDge);
			for (int i = 0; i < EDge; i++)
				capacity[i] = G.incL[i].capacity;
			chormes = (int**)malloc(sizeof(int*)*pop);
			childs = (int**)malloc(sizeof(int*)*Beta);
			monsters = (int**)malloc(sizeof(int*)*Gama);
			for (int i = 0; i < pop; i++)
				chormes[i] = (int*)malloc(sizeof(int)*Task);
			for (int i = 0; i < Beta; i++)
				childs[i] = (int*)malloc(sizeof(int)*Task);
			for (int i = 0; i<Gama; i++)
				monsters[i] = (int*)malloc(sizeof(int)*Task);
		};
		~NewGA(){
			free(st);
			free(te);
			free(demand);
			free(capacity);
			free(chormes);
			free(childs);
			free(monsters);
		}
	public:

		//vector<int>
		vector<int> GAsearch(vector<service>&nser, taskPath*PathSets){
			cout<<"ops"<<endl;
			std::vector<service>ser;
			int c=0;
			std::vector<int> flag(nser.size(),-1);
			for(int i=0;i<nser.size();i++)
			{	
				int co=0;
				if(PathSets[i].num>=4)
					while(PathSets[i].Pathset[co][0]>=0&&co<4)co++;
				if(i!=77&&i!=168&&co>=4)
					{	
						flag[i]=1;
						ser.push_back(nser[i]),c++;
					}
			}
			int T=c;
			int  W=4;
			taskPath*_PathSets=new taskPath[T];
			int p=0;
			for(int i=0;i<nser.size();i++)
					if(flag[i]>0)
							_PathSets[p++]=PathSets[i];

			vector<double>capacity(EDge,100);
			vector<pair<string,double>> rdata;
			vector<double>u(EDge,0.1);
			vector<vector<vector<double>>>spu(ser.size(),std::vector<std::vector<double>>());
			vector<double>f(EDge,0);
			vector<int>ids(ser.size(),0);
			vector<int>bids(ser.size(),0);
			vector<vector<pair<int,int>>>epath(EDge,vector<pair<int,int>>());
			map<pair<int,int>,vector<int>>pathe;
			vector<vector<double>>x(ser.size(),vector<double>());
			vector<vector<double>>cc(ser.size(),vector<double>());
			vector<double>y(ser.size(),0);
			vector<vector<double>>spy(ser.size(),vector<double>());
			vector<vector<double>>overtime(ser.size(),std::vector<double>());
			map<pair<int,int>,double>wsigu;
			map<pair<int,pair<int,int>>,double>overo;
			vector<double> tsum(ser.size(),-1e10);
			vector<double>esigxt(EDge,0);
			double alpha=6;
			double q=1.0;//0.98;
			std::vector<double>keep(ser.size(),1000);
			cout<<"asd"<<endl;
			cout<<"service size is "<<ser.size()<<endl;
			for(int i=0;i<ser.size();i++)
			{
				double as=0;
				int num=min(_PathSets[i].num,4);
				int j=0;
				if(_PathSets[i].Pathset[j][0]>-1)
					{
						std::vector<double> v;
						x[i].push_back(0);
						spy[i].push_back(100);
						spu[i].push_back(v);
						for(int M=0;M<EDge;M++)
							spu[i][spu[i].size()-1].push_back(0.1);
						as+=100;
					}
				int k=0;
				while(_PathSets[i].Pathset[j][k]>-1)
					{
						epath[_PathSets[i].Pathset[j][k]].push_back(make_pair(i,0));
						pathe[make_pair(i,0)].push_back(_PathSets[i].Pathset[j][k]);
						k++;
					}
					cc[i].push_back(0);
					overtime[i].push_back(1);
				for(int jj=1;jj<num;jj++)
					{
						int j=jj;
						if(_PathSets[i].Pathset[j][0]>-1)
							{
								std::vector<double> v;
								x[i].push_back(0);
								spy[i].push_back(100);
								spu[i].push_back(v);
								for(int M=0;M<EDge;M++)
									spu[i][spu[i].size()-1].push_back(0.1);
								as+=100;
							}
						int k=0;
						while(_PathSets[i].Pathset[j][k]>-1)
							{
								epath[_PathSets[i].Pathset[j][k]].push_back(make_pair(i,jj));
								pathe[make_pair(i,jj)].push_back(_PathSets[i].Pathset[j][k]);
								k++;
							}
							cc[i].push_back(0);
							overtime[i].push_back(1);
					}
				y[i]=100;
			}

			for(int i=0;i<ser.size();i++)
				for(int j=0;j<x[i].size();j++)
					if(pathe[make_pair(i,j)].size()<=0)
						cout<<"fucker"<<endl;
			double bestv=-INT_MAX/2;
			int count=0;
			double tt=1;
			double beta=0.05;
			double gama=0.05;
			int counter=0;
			time_t begin=clock();
			for(int t=0;t<100000;t++)
			{
				vector<double> sigxt(ser.size(),0);
				map<pair<int,int>,double>sigu;
				vector<double>asigxt(ser.size(),0);
				double objective=0;
				for(int i=0;i<EDge;i++)
					{
							if(esigxt[i]-capacity[i]>1)
								{
									counter++;
									if(counter>100000&&beta>0){
										//beta=0;
										//gama/=5;
										//cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
									}
									u[i]=u[i]+gama*u[i]*(esigxt[i]-capacity[i])/capacity[i];
							}	
							else
								u[i]=u[i]+beta*u[i]*(esigxt[i]-capacity[i])/capacity[i];
					}
				for(int i=0;i<EDge;i++)
					{
						esigxt[i]=0;
						for(int j=0;j<epath[i].size();j++)
							sigu[epath[i][j]]+=u[i];

					}
				for(int i=0;i<ser.size();i++)
				{
					//for(int j=0;j<x[i].size();j++)
					asigxt[i]=x[i][ids[i]];
					y[i]=y[i]+(0.05/12)*(asigxt[i]-y[i]);
				}
				for(int i=0;i<ser.size();i++)
					{
						double min=DBL_MAX;
						int id=-1;
						double toto=0;
						double sum=0;
						for(int j=0;j<x[i].size();j++)
						{
							double data=sigu[make_pair(i,j)];
							if(data<min){
								min=data;
								id=j;
							}
						}
						for(int j=0;j<x[i].size();j++)
							x[i][j]=0;
						x[i][id]=pow(pow(y[i],-alpha)/sigu[make_pair(i,id)],10)*y[i];
						for(int k=0;k<pathe[make_pair(i,id)].size();k++)
							esigxt[pathe[make_pair(i,id)][k]]+=x[i][id];
						ids[i]=id;
					}
				int overflow=0;
				for(int i=0;i<ser.size();i++)
					objective+=-pow(x[i][ids[i]],-5)/5;
				if(objective>bestv&&overflow<100)
				{
					bestv=objective;
					bids=ids;
				}
				for(int k=0;k<EDge;k++)
				{
					//cout<<esigxt[k]<<" ";
			//kk+=pow(x[k],-5)/5;
				}
				//cout<<endl;
				//for(int k=0;k<ser.size();k++)
					//cout<<pow(x[k][ids[k]],-5)/5<<" ";
				//cout<<endl;
				
				cout<<"obj is: "<<t<<" "<<objective<<" "<<counter<<endl;
			}
			time_t end=clock();
			cout<<"time is "<<end-begin<<endl;
			return bids;
		};
	};