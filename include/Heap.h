#include<vector>
#define N 100000
using namespace std; 
struct edge
{
public:
	int head, tail;
    double weight;
	float capacity;
	float backweight;
	edge(int _head,double _weight,float _capacity=0,float _backweight=0):head(_head),tail(0),weight(_weight),capacity(_capacity),backweight(_backweight){
	};
};
class Heap
{
public:
	Heap(int n);
	~Heap();
	void push(int vertID, double w);
	void update(int vertID, double w);
	int pop();
	int empty();
private:
	vector<edge>h;
	vector<int> post;
	int nodeNum;
	void fix(int fixID);
};


