#include <set>
#include <vector>
#include <string>
#include <utility>
#include <iostream>

#include <random>
#include <chrono>

#include <cstdio>
#include <cstdlib>

#include "network.hpp"
#include "diffusionstate.hpp"
#include "utils.hpp"
#include "tools.hpp"

using namespace std;

void testUpLow(const Network &network, DiffusionState_MIC &diffu, vector<rTuple> &rtup, int *nodes){
    set<int> seedset;
    for(int i = 0;i < 4;i++){
        set<int> seed;
        for(int j = i*100;j < i*100+100;i++)
            seed.insert(nodes[i]);
        diffu.seed(seed);
    }
    for(int i = 1000;i < network.vertexNum;i += 5000)
        seedset.insert(nodes[i]);
    
    int cindex = diffu.seed(seedset);
    cout << diffu.expInfluenceComplete(network, 5000, cindex) << endl;
    diffu.removeSeed(cindex);
    cout << diffu.computeG(seedset, rtup, network.vertexNum, "upper", nullptr) << endl;
    cout << diffu.computeG(seedset, rtup, network.vertexNum, "mid", nullptr) << endl;
    cout << diffu.computeG(seedset, rtup, network.vertexNum, "lower", nullptr) << endl;
}

int main(int args, char **argv){
    string name = string(argv[1]);
    string type = string(argv[2]);
    int vnum = atoi(argv[3]);
    string priority = string(argv[4]);

    string path = "../data/"+name+".txt", fname;

    Network network(path, type, vnum);
    network.setICProb(0.1);

    mt19937 rand(chrono::high_resolution_clock::now().time_since_epoch().count());
    DiffusionState_MIC diffusionState(network, priority, rand);

    int n, shuffle_node[vnum];
    path = "../data/"+name+"_node.txt";
    FILE *fd = fopen(path.c_str(), "r");
    for(int i = 0;i < vnum;i++){
        fscanf(fd, "%d", &n);
        shuffle_node[i] = n;
    }
    fclose(fd);

    vector<rTuple> rtup, all_rtup;
    set<int> check;
    path = "../data/"+name+"_rtup.txt";
    fd = fopen(path.c_str(), "r");
    for(int i = 0;i < 5000000;i++){
        rTuple rt;
        rt.retrieve(fd);
        all_rtup.push_back(rt);
    }
    fclose(fd);
    int pick;
    while(check.size() < 2000000){
        pick = rand() % 5000000;
        if(check.find(pick) != check.end())
            continue;
        rtup.push_back(all_rtup[pick]);
    }

    testUpLow(network, diffusionState, rtup, shuffle_node);

    return 0;
}