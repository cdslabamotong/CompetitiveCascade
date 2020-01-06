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

#define NUM_RTUP 5000000

using namespace std;

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

    vector<rTuple> rtup;
    path = "../data/"+name+"_rtup.txt";
    FILE *fd = fopen(path.c_str(), "w");

    diffusionState.getRTuples(network, rtup, NUM_RTUP);

    for(int i = 0;i < NUM_RTUP;i++)
        rtup[i].store(fd);
    
    fclose(fd);
    
    return 0;
}