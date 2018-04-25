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
__device__ double Add(double* address, double val)
{
    unsigned long long int* address_as_ull =
                             (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
__global__ void getpath(Edge*edge,int*p,double*x,double*f,int*s,int*t,int NN,int T)
{
	int id=threadIdx.x + blockDim.x*blockIdx.x;
	if(id>=T)return;
	int pre=tid[id];
	int off=NN*s[id];
	int ss=s[id];
	double v=0;
	int tm=-1;
	while(pre!=ss)
		{
			tm=p[pre];
			pre=edge[tm].head;
			Add(&f[tm],x[id]);
			
			
		}
	
	
	
	
}
__global__ void ChangePameterC(int*p,double*d,int n){
	int tid = blockIdx.y;
	int i = threadIdx.x + blockDim.x*blockIdx.x;
	if (i>=n||tid >=n)return;
	int biao = tid*n + i;
	d[biao] = (i == tid) ? 0.0:DBL_MAX/2;
	p[biao] = -1;
}
__global__ void bellmanHigh(Ldge*edge, int*m, double*c, int*p, double*u,int E,int NN)
{
	int tid = blockIdx.y;
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if (i >=E)return;
	int head = edge[i].head;
	int tail = edge[i].tail;
	int biao = tid*NN+head;
	float val = c[tid*NN+tail]+u[i];
	if (c[biao] >val){
		*m = 1;
		c[biao] = val;
	}
}
__global__ void color(Ldge *edge, int *m, double*c,int*p,double*u,int E,int NN){

	int tid = blockIdx.y;
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if (i >=E)return;
	int head = edge[i].head;
	int tail = edge[i].tail;
	int biao = tid*NN+head;
	float val = c[tail+tid*NN]+1+u[i];
	if (c[biao] == val){
		p[biao] = tid;
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
	cudaMalloc((void**)&dev_sum, (T/512+1)*sizeof(double));
	cudaMalloc((void**)&dev_d, NN*NN*sizeof(double));
	cudaMalloc((void**)&dev_p, NN*NN*sizeof(int));
	cudaMalloc((void**)&dev_edge, E*sizeof(Ldge));
	cudaMemcpy(dev_x,x, sizeof(double)*(T+1),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_y,y, sizeof(double)*T,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_u,u, sizeof(double)*E,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_f,f, sizeof(double)*E,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_sum,sum,sizeof(double)*(T/512+1),cudaMemcpyHostToDevice);
	cudaMemcpy(dev_d,d,sizeof(double)*NN*NN,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_p,p,sizeof(int)*NN*NN,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_edge,edge,sizeof(Ldge)*E,cudaMemcpyHostToDevice);
}
vector<pair<string,float> > NewGAParrel::GAsearch(){
	Cudamalloc();
	cout<<"ajsgkvksafuaqd"<<endl;
	time_t begin=clock();
	for(int i=0;i<100000;i++)
	{
		dim3 blocksq(NN/64+1, NN*NN);
		ChangePameterC << <blocksq,64>> >(dev_p, dev_d,NN);
		dim3 blocks_square(E/256+1,NN);
		do{
			cudaMemcpy(dev_m,m, sizeof(int), cudaMemcpyHostToDevice);
			bellmanHigh << <blocks_square,256>> >(dev_edge, dev_m, dev_d, dev_p, dev_u,E,NN);
			cudaMemcpy(m, dev_m, sizeof(int), cudaMemcpyDeviceToHost);
		} while (*m);
		
	}
	time_t end=clock();
	cout<<end-begin<<endl;
	cout<<"what happened"<<endl;
}
