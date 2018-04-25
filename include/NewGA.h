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
#include"Heap.h"
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
	//std::vector<map<int,int>> serid;
	std::vector<vector<pair<int,int>>>gugu;
	std::vector<set<int>> sset;
	public:
		NewGA( Graph &_G):G(_G){
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
			cout<<"hjk: "<<nser.size()<<endl;
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
									if(counter>200000&&beta>0){
										beta=0;
										gama/=5;
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
					asigxt[i]=0;
					for(int j=0;j<x[i].size();j++)
						asigxt[i]+=x[i][j];
					y[i]=y[i]+(0.1/12)*(asigxt[i]-y[i]);
				}
				double ttf=0;
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
							x[i][j]=pow(pow(y[i],-alpha)/sigu[make_pair(i,j)],10)*y[i];
							for(int k=0;k<pathe[make_pair(i,j)].size();k++)
								esigxt[pathe[make_pair(i,j)][k]]+=x[i][j];
							ttf+=x[i][j];
						}
						/*for(int j=0;j<x[i].size();j++)
							x[i][j]=0;
						x[i][id]=pow(pow(y[i],-alpha)/sigu[make_pair(i,id)],10)*y[i];
						ttf+=x[i][id];
						for(int k=0;k<pathe[make_pair(i,id)].size();k++)
							esigxt[pathe[make_pair(i,id)][k]]+=x[i][id];
						ids[i]=id;*/
					}
				cout<<"total flow "<<ttf<<endl;
				int overflow=0;
				for(int i=0;i<ser.size();i++)
					objective+=-pow(asigxt[i],-5)/5;
				if(objective>bestv&&overflow<100)
				{
					bestv=objective;
					bids=ids;
				}
				cout<<"obj is: "<<t<<" "<<objective<<" "<<counter<<endl;
			}
			time_t end=clock();
			cout<<"time is "<<end-begin<<endl;
			return bids;
		};
		pair<vector<int>,double> diks(int s,int t,std::vector<double>&v)
		{
			int n_num = G.n;
			vector<int> flag(n_num,0);
			vector<int> peg(n_num,-1);
			vector<double>d(n_num,0);
			for (int i = 0; i < n_num; i++)
				if (i == s)
					d[i] = 0.0;
				else
					d[i] =DBL_MAX/2;
			Heap heap(n_num);
			for (int i = 0; i < n_num; i++)
				heap.push(i, d[i]);
			do{
				int cur=heap.pop();
				if(cur==t)break;
				if(flag[cur]==1)continue;
				flag[cur] = 1;
				int size = G.near[cur].size();
				//cout<<cur<<endl;
				for (int i = 0; i<size; i++){
					int id = G.near[cur][i];	
					Edge* e = &G.incL[id];
					int delt = 0;
					if (d[e->head]>d[e->tail]+v[id]){
						//cout<<d[e->head]<<" ";
						d[e->head] = d[e->tail] + v[id];
					//	cout<<d[e->head]<<" "<<d[e->tail]<<" "<<v[id]<<endl;
						heap.update(e->head, d[e->head]);
						peg[e->head] = id;
					}
				}
			} while (!heap.empty());


			//cout<<"out "<<endl;
			vector<int> rout;
			int p=t;
			double value=0;
			//cout<<"outmmmmmmmmmmmmmmmm "<<endl;
			while(p!=s)
			{
				rout.push_back(peg[p]);
				value+=v[peg[p]];
				p=G.incL[peg[p]].tail;
			}
			return make_pair(rout,value);
		};
		void diksset(vector<vector<double>>&x,vector<double>&f,int s,set<int>&sset,std::vector<double>&v,vector<pair<int,int>> ar)
		{
			int sizee=ar.size();
			int n_num = G.n;
			vector<int> flag(n_num,0);
			vector<int> peg(n_num,-1);
			vector<double>d(n_num,0);
			for (int i = 0; i < n_num; i++)
				if (i == s)
					d[i] = 0.0;
				else
					d[i] =DBL_MAX/2;
			Heap heap(n_num);
			for (int i = 0; i < n_num; i++)
				heap.push(i, d[i]);

			do{
				int cur=heap.pop();
				if(flag[cur]==1)continue;
				if(sset.find(cur)!=sset.end()&&flag[cur]!=1)sizee--;
				flag[cur]=1;
				if(sizee==0)break;
				flag[cur] = 1;
				int size = G.near[cur].size();
				for (int i = 0; i<size; i++){
					int id = G.near[cur][i];	
					Edge* e = &G.incL[id];
					int delt = 0;
					if (d[e->head]>d[e->tail]+v[id]){
						d[e->head] = d[e->tail] + v[id];
						heap.update(e->head, d[e->head]);
						peg[e->head] = id;
					}
				}
			} while (!heap.empty());
			double value=0;
			vector<pair<vector<int>,double>>re;
			int iter;

			for(int iter=0;iter<ar.size();iter++)
			{
				double value=0;
				int p=ar[iter].first;
				while(p!=s)
				{
					value+=v[peg[p]];
					f[peg[p]]+=x[ar[iter].second][0];
					p=G.incL[peg[p]].tail;
				}
				x[ar[iter].second][0]=pow(value,-(1.0/6.0));
				
			}
		};
	void getpath(vector<vector<double>>&x,vector<double>&u,std::vector<double>&flow)
	{
		for(int s=0;s<sset.size();s++)
			diksset(x,flow,s,sset[s],u,gugu[s]);

	};
	vector<int> GAsearchP(vector<service>&nser, taskPath*PathSets){
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
			std::vector<set<int>> tmps(G.n,set<int>());
			sset=tmps;
			std::vector<std::map<int, int>> vv;
			std::vector<vector<pair<int,int>>> gu(G.n,vector<pair<int,int>>());
			gugu=gu;
			for(int i=0;i<ser.size();i++)
			{
				sset[ser[i].s].insert(ser[i].t);
				gugu[ser[i].s].push_back(make_pair(ser[i].t,i));

			}
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
			vector<vector<int>>serpath(ser.size(),vector<int>());
			vector<double>esigxt(EDge,0);

			vector<double> tsum(ser.size(),-1e10);
			double alpha=6;
			double q=1;
			std::vector<double>keep(ser.size(),1000);
			cout<<"asd"<<endl;
			for(int i=0;i<ser.size();i++)
			{
				double as=0;
				int num=min(_PathSets[i].num,1);
				int j=rand()%num;
				if(_PathSets[i].Pathset[j][0]>-1)
					{
						std::vector<double> v;
						x[i].push_back(100);
						spy[i].push_back(100);
						spu[i].push_back(v);
						for(int M=0;M<EDge;M++)
							spu[i][spu[i].size()-1].push_back(0.1);
						as+=100;
					}
				int k=0;
				vector<int> rout;
				while(_PathSets[i].Pathset[j][k]>-1)
					{
						rout.push_back(_PathSets[i].Pathset[j][k]);
						epath[_PathSets[i].Pathset[j][k]].push_back(make_pair(i,0));
						pathe[make_pair(i,0)].push_back(_PathSets[i].Pathset[j][k]);
						k++;
					}
					serpath[i]=rout;
					cc[i].push_back(0);
					overtime[i].push_back(1);
				for(int jj=1;jj<num;jj++)
					{
						int j=rand()%num;
						if(_PathSets[i].Pathset[j][0]>-1)
							{
								std::vector<double> v;
								x[i].push_back(100);
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
			time_t start=clock();
			time_t toto=0;
			for(int t=0;t<100000;t++)
			{
				int over=0;
				for(int i=0;i<EDge;i++)
					{
							if(esigxt[i]-capacity[i]>1)
								{
									over++;
									counter++;
									if(counter>200000&&beta>0){
										beta=0;
										gama/=10;
									}
									u[i]=max(u[i]+gama*u[i]*(esigxt[i]-capacity[i])/capacity[i],0.0);
							}	
							else
								u[i]=max(u[i]+beta*u[i]*(esigxt[i]-capacity[i])/capacity[i],0.0);
							esigxt[i]=0;
					}
				time_t s=clock();
				getpath(x,u,esigxt);
				time_t e=clock();
				toto+=(e-s);
				double objective=0;
				for(int i=0;i<ser.size();i++)
					objective+=-pow(x[i][0],-5)/5;
				cout<<"obj is: "<<over<<" "<<" "<<objective<<" "<<counter<<endl;
			}
			time_t end=clock();
			cout<<toto<<" "<<end-start<<endl;
			return bids;
		};
	};