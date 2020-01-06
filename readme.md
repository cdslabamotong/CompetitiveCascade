# On Multi-Cascade Influence Maximization: Model, Hardness and Algorithmic Framework

This page shows the accompany C++ code for online paper _[On Multi-Cascade Influence Maximization: Model, Hardness and Algorithmic Framework](https://arxiv.org/abs/1912.00272)_

## Data Structure

The code entry locates in the file `start.cpp`, and the dataset will be read and stored in the class `Network` in `network.hpp`.

```cpp
class Network{
  typedef struct _Edge {
    int u, v;
    double p;
  } _Edge;    // structure for storing edges
  
  ...

public:
  int vertexNum, edgeNum;
  int *neighbor, *neighbor_reverse;
  double *probability;

  vector<int> sorted_degree;
  int *neighbor_ptr, *neighbor_reverse_ptr;
  int *inDegree, *outDegree;
  int *s_contri_order;
  double *threshold, *c_threshold, *s_contri;
  string type, path, s_contri_path;
  double IC_prob;
  bool is_s_contri;
  
  ...
  
};
```

During the process of diffusion, the class `DiffusionState_MIC` in `diffusionstate.hpp` maintains the status of nodes and graph. Besides, class `DiffusionState_MIC` will control the diffusion policy and priority.

```cpp
class DiffusionState_MIC{
  std::set<int> new_active[THREAD];

  ...

  int priority(const std::set<std::pair<int,short>> &pri, int n, int tid, int shift=0){
    switch(type_bak[0]){
    case 'c':
      return pri_cas(pri, n, shift);
    case 'r':
      return pri_rand(pri, tid);
    case 'n':
      return pri_nei(pri);
    }
  }

  ...

public:
  short *seed_state;
  std::set<int> seednodes;
  std::map<int, std::set<int> > seedsets;
  
  ...

  void diffuse(const Network &network, int *result, int cindex, int round, int j){
    
    ...
    
    new_active[tid] = seednodes;
    int base = vnum * tid;
    for(int i = 0;i < vnum;i++)
      temp_state_1[base+i] = seed_state[i];
    for(int i = 0;i < round;i++){
      DiffusionState_MIC::diffuseOneRound(network, temp_state_1+base, tid);
      if(new_active[tid].empty())  break;
    }
    for(int i = 0;i < vnum;i++)
      if(temp_state_1[base+i] == cindex)
        result[j]++;
    
    ...
    
  }
  
  ...
  
};
```

## Requirement

Our code works based on _Boost_ with version later than _1.68.0_

## Compile & Run

To compile the code, first make sure _Boost_ library is installed correctly, and then use following command

```bash
$ g++ -o start start.cpp -lboost_thread -lboost_system -lpthread -std=c++11
```

Or, for convenience

```bash
$ make
```

After the executable file is created, it could be run using command

```bash
$ ./start <dataset> <diffusion policy: WC/IC/VIC> <num of nodes> <num of seed k> <test span> <diffusion priority: random/cascade/neighbor> <load result? (default: false)> <new test? (default: true)>
```
