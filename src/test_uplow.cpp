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
    set<int> seedset, check;
    int n;
    while(check.size() < 100){
        n = rand() % (network.vertexNum-400);
        if(check.find(n) != check.end())
            continue;
        check.insert(n);
        seedset.insert(nodes[n+400]);
    }
    cout << seedset.size() << endl;
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

    string path = "../data/"+name+".txt";
    char fname[256];

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

    vector<rTuple> rtup;
    chrono::high_resolution_clock::time_point start, end;
    set<int> seedset;
    for(int i = 0;i < 4;i++){
        seedset.clear();
        for(int j = 100*i;j < 100*i+100;j++)
            seedset.insert(shuffle_node[j]);
        diffusionState.seed(seedset);
    }
    for(int i = 0;i < 50;i++){
        rtup.clear();
        start = chrono::high_resolution_clock::now();
        diffusionState.getRTuples(network, rtup, 2000000);
        testUpLow(network, diffusionState, rtup, shuffle_node);
        end = chrono::high_resolution_clock::now();
        printTime(start, end);
        cout << "----------" << endl;
    }

    return 0;
}