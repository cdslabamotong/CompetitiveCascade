#pragma once

#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <algorithm>

#include "tools.hpp"

class rTuple{
public:
    int node_v;
    bool isdiff;
    std::set<int> upper, lower, seed;
    std::map<int,std::set<int>> relations;
    rTuple(){clear();}
    rTuple(const rTuple &rt){
        node_v = rt.node_v;
        isdiff = rt.isdiff;
        upper = rt.upper;
        lower = rt.lower;
        relations = rt.relations;
    }
    rTuple &operator = (const rTuple &)=delete;

    void clear(){
        node_v = -1;
        upper.clear();
        lower.clear();
        seed.clear();
        relations.clear();
        isdiff = false;
    }

    void _stat(){
        std::cout << "----------" << std::endl;
        std::cout << "node_v: " << node_v << std::endl;
        std::cout << "upper:" << std::endl;
        printContainer(upper);
        std::cout << "lower:" << std::endl;
        printContainer(lower);
        std::cout << "----------" << std::endl;
    }

    void store(FILE *fd){
        fprintf(fd, "%d %d", node_v, isdiff?1:0);
        fprintf(fd, "\n%ld", upper.size());
        for(int i: upper)
            fprintf(fd, " %d", i);
        fprintf(fd, "\n%ld", lower.size());
        for(int i: lower)
            fprintf(fd, " %d", i);
        fprintf(fd, "\n%ld", seed.size());
        for(int i: seed)
            fprintf(fd, " %d", i);
        fprintf(fd, "\n%ld", relations.size());
        for(std::pair<int,std::set<int>> i: relations){
            fprintf(fd, "\n%d\n%ld ", i.first, i.second.size());
            for(int j: i.second)
                fprintf(fd, " %d", j);
        }
        fprintf(fd, "\n");
    }

    void retrieve(FILE *fd){
        upper.clear();
        lower.clear();
        seed.clear();
        relations.clear();
        int cache, count;
        fscanf(fd, "%d %d", &node_v, &cache);
        isdiff = cache?true:false;
        fscanf(fd, "%d", &count);
        for(int i = 0;i < count;i++){
            fscanf(fd, "%d", &cache);
            upper.insert(cache);
        }
        fscanf(fd, "%d", &count);
        for(int i = 0;i < count;i++){
            fscanf(fd, "%d", &cache);
            lower.insert(cache);
        }
        fscanf(fd, "%d", &count);
        for(int i = 0;i < count;i++){
            fscanf(fd, "%d", &cache);
            seed.insert(cache);
        }
        int ssize, idx;
        fscanf(fd, "%d", &count);
        for(int i = 0;i < count;i++){
            fscanf(fd, "%d", &idx);
            fscanf(fd, "%d", &ssize);
            relations[idx] = std::set<int>();
            for(int j = 0;j < ssize;j++){
                fscanf(fd, "%d", &cache);
                relations[idx].insert(cache);
            }
        }
    }
};