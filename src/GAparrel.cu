#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include"Graph.h"
#include "service.h"
#include"taskPath.h"
#include"valuemark.h"
#include"curand_kernel.h"
#include"iostream"
#include <fstream>
#include"const.h"
#include<math.h>
#include"BFS.h"
#include"GAparrel.h"
#include<time.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/functional.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include"PathArrange.h"
__global__ void PathChoose(int T,int M,int W,double *x,double *y,double*u,double*f,int*paths)
{
	int id = threadIdx.x + blockIdx.x*blockDim.x;
	int tid=threadIdx.x;
	if(id>=T*W)return;
	int d=id/W;
	int k=id%W;
	__shared__ double price[256];
	price[tid]=0;
	int off=d*M*W+k*M;
	for(int i=0;i<M;i++)
		{
		if(paths[off+i]<0)break;
		price[tid]+=u[paths[off+i]];
		}
	if(price[tid]<=0)price[tid]=DBL_MAX;
	if(k==0)
	{
		double gu=price[tid];
		int bid=0;
		for(int i=1;i<W;i++)
			if(price[tid+i]<gu)
				{
				gu=price[tid+i];
				bid=i;
				}
		y[d]+=(0.05/12)*(x[d]-y[d]);
		x[d]=pow(pow(y[d],-6)/price[tid+bid],10)*y[d];
		int off=d*M*W+bid*M;
		for(int i=0;i<M;i++)
			{
			if(paths[off+i]<0)break;
			f[paths[off+i]]+=x[d];
			}
	}
}
__global__ void changeU(int E,double*u,double*f)
{
	int id = threadIdx.x + blockIdx.x*blockDim.x;
	if(id>=E)return;
	u[id]+=u[id]*0.05*(f[id]-100.0)/100.0;
	f[id]=0;
}
__global__ void Sum(int T,double*x,double*sum)
{
	int id = threadIdx.x + blockIdx.x*blockDim.x;
	if(id>T/2)return;
	int bid=blockIdx.x;
	int off=bid*1024;
	int left=min(T-512*bid,512);
	int tid=threadIdx.x;
	__shared__ double tsum[256];
	tsum[tid]=pow(x[tid+off],-5)/5;
	if(id>=T/2)return;
	tsum[tid]+=pow(x[tid+off+(left+1)/2],-5)/5;
	for(int s=(left+1)/2;s>1;s=(s+1)/2)
	{
		if(tid<s/2)
			tsum[tid]+=tsum[tid+(s+1)/2];
		__syncthreads();
	}
	__syncthreads();
	if(tid==0)
		sum[bid]=tsum[0];
		
}
void NewGAParrel::Cudamalloc(){
	cout<<"m is"<<M<<endl;
	cudaMalloc((void**)&dev_x, (T+1)*sizeof(double));
	cudaMalloc((void**)&dev_y, T*sizeof(double));
	cudaMalloc((void**)&dev_u, E*sizeof(double));
	cudaMalloc((void**)&dev_f, E*sizeof(double));
	cudaMalloc((void**)&dev_paths, T*M*W*sizeof(int));
	cudaMalloc((void**)&dev_sum, (T/512+1)*sizeof(double));
	cudaMemcpy(dev_x,x, sizeof(double)*(T+1),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y,y, sizeof(double)*T,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_u,u, sizeof(double)*E,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_f,f, sizeof(double)*E,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_paths,paths,sizeof(int)*T*W*M,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_sum,sum,sizeof(double)*(T/512+1),cudaMemcpyHostToDevice);
}
vector<pair<string,float> > NewGAParrel::GAsearch(){
	Cudamalloc();
	cout<<"m is "<<M<<endl;
	for(int i=0;i<100000;i++)
	{
		PathChoose<< <T*W/256+1,256>> >(T,M,W,dev_x,dev_y,dev_u,dev_f,dev_paths);
		cudaMemcpy(f,dev_f, sizeof(double)*E,cudaMemcpyDeviceToHost);
		changeU<< <E/512+1,512>> >(E,dev_u,dev_f);
		Sum<<<T/256+1,256>>>(T,dev_x,dev_sum);
		cudaMemcpy(sum,dev_sum,sizeof(double)*(T/1024+1),cudaMemcpyDeviceToHost);
		cout<<"sum o is: "<<sum[0]<<endl;
	}
	cout<<"what happened"<<endl;
}
