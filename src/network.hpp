#pragma once

#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <map>
#include <set>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

using namespace std;

class Network{
    typedef struct _Edge {
        int u, v;
        double p;
    } _Edge;

    void clearPointers(){
        neighbor = nullptr;
        neighbor_reverse = nullptr;
        neighbor_ptr = nullptr;
        neighbor_reverse_ptr = nullptr;
        probability = nullptr;
        inDegree = nullptr;
        outDegree = nullptr;
        s_contri_order = nullptr;
        threshold = nullptr;
        c_threshold = nullptr;
        s_contri = nullptr;
    }

    template <typename T>
    void freeSpace(T *ptr){
        if(ptr != nullptr)
            delete [] ptr;
    }

    void importRelation(string path){
        _Edge *edge = new _Edge[100000000];
        char line[50];
        int u, v, tmp;
        // set<int> nodeSet;
        
        FILE *fd = fopen(path.c_str(), "r");
        while(NULL != fgets(line, 50, fd)){
            vector<string> inStr;
            boost::split(inStr, line, boost::is_any_of(" "));
            u = stoi(inStr[0]);
            v = stoi(inStr[1]);
            edge[edgeNum].u = u;
            edge[edgeNum].v = v;
            inDegree[v]++;
            outDegree[u]++;
            if(type == "VIC")
                edge[edgeNum].p = stod(inStr[2]);
            edgeNum++;
        }
        fclose(fd);

        neighbor = new int[edgeNum];
        neighbor_reverse = new int[edgeNum];
        probability = new double[edgeNum];
        sort(edge, edge+edgeNum, [](const _Edge &a, const _Edge &b)-> bool {
            if(a.u == b.u)
                return a.v < b.v;
            return a.u < b.u;
        });
        for(int i = 0;i < edgeNum;i++){
            neighbor[i] = edge[i].v;
            probability[i] = edge[i].p;
        }

        sort(edge, edge+edgeNum, [](const _Edge &a, const _Edge &b)-> bool {
            if(a.v == b.v)
                return a.u < b.u;
            return a.v < b.v;
        });
        for(int i = 0;i < edgeNum;i++)
            neighbor_reverse[i] = edge[i].u;

        memcpy(neighbor_ptr, outDegree, vertexNum*sizeof(int));
        memcpy(neighbor_reverse_ptr, inDegree, vertexNum*sizeof(int));
        for(int i = 1;i < vertexNum;i++){
            neighbor_ptr[i] += neighbor_ptr[i-1];
            neighbor_reverse_ptr[i] += neighbor_reverse_ptr[i-1];
        }
        
        delete [] edge;
    }

    inline double getProbability(int i, int j) const {
        if(i)   return probability[neighbor_ptr[i-1]+j];
        return probability[j];
    }

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

    Network(){
        clearPointers();
        is_s_contri = false;
    }

    Network(string path, string type, int vertexNum){
        clearPointers();
        cout << "import " << path << " " << type << endl;
        this->path = path;
        this->type = type;
        this->vertexNum = vertexNum;
        this->edgeNum = 0;
        is_s_contri = false;
        inDegree = new int[vertexNum];
        outDegree = new int[vertexNum];
        neighbor_ptr = new int[vertexNum];
        neighbor_reverse_ptr = new int[vertexNum];
        sorted_degree.clear();

        importRelation(path);

        if(type == "IC")
            IC_prob = .1;
        else if(type == "LT"){
            threshold = new double[vertexNum];
            c_threshold = new double[vertexNum];
        }
        else if(type != "WC" && type != "VIC")
            printf("Invalid model\n");
    }

    Network(const Network &network){
        clearPointers();
        path = network.path;
        type = network.type;
        vertexNum = network.vertexNum;
        edgeNum = network.edgeNum;
        is_s_contri = network.is_s_contri;
        neighbor = new int[edgeNum];
        neighbor_reverse = new int[edgeNum];
        memcpy(neighbor, network.neighbor, edgeNum*sizeof(int));
        memcpy(neighbor_reverse, network.neighbor_reverse, edgeNum*sizeof(int));
        neighbor_ptr = new int[vertexNum];
        neighbor_reverse_ptr = new int[vertexNum];
        memcpy(neighbor_ptr, network.neighbor_ptr, vertexNum*sizeof(int));
        memcpy(neighbor_reverse_ptr, network.neighbor_reverse_ptr, vertexNum*sizeof(int));
        probability = new double[edgeNum];
        memcpy(probability, network.probability, edgeNum*sizeof(double));
        sorted_degree = network.sorted_degree;
        IC_prob = network.IC_prob;
        s_contri_path = network.s_contri_path;
        inDegree = new int[vertexNum];
        outDegree = new int[vertexNum];
        memcpy(inDegree, network.inDegree, vertexNum*sizeof(int));
        memcpy(outDegree, network.outDegree, vertexNum*sizeof(int));
        if(network.s_contri != nullptr){
            s_contri_order = new int[vertexNum];
            s_contri = new double[vertexNum];
            memcpy(s_contri_order, network.s_contri_order, vertexNum*sizeof(int));
            memcpy(s_contri, network.s_contri, vertexNum*sizeof(double));
        }
        if(type == "LT"){
            threshold = new double[vertexNum];
            c_threshold = new double[vertexNum];
            memcpy(threshold, network.threshold, vertexNum*sizeof(double));
            memcpy(c_threshold, network.c_threshold, vertexNum*sizeof(double));
        }
    }

    ~Network(){
        freeSpace(neighbor);
        freeSpace(neighbor_reverse);
        freeSpace(neighbor_ptr);
        freeSpace(neighbor_reverse_ptr);
        freeSpace(probability);
        freeSpace(inDegree);
        freeSpace(outDegree);
        freeSpace(s_contri);
        freeSpace(s_contri_order);
        freeSpace(threshold);
        freeSpace(c_threshold);
    }

    inline int getNeighbor(int i, int j) const {
        if(i)   return neighbor[neighbor_ptr[i-1]+j];
        return neighbor[j];
    }

    inline int getReverseNeighbor(int i, int j) const {
        if(i)   return neighbor_reverse[neighbor_reverse_ptr[i-1]+j];
        return neighbor_reverse[j];
    }

    void setICProb(double prob){
        IC_prob = prob;
    }

    void sortByDegree(){
        vector<pair<int,int>> tempMap;
        for(int i = 0;i < vertexNum;i++)
            tempMap.push_back(make_pair(i,outDegree[i]));
        sort(tempMap.begin(), tempMap.end(), [](const pair<int,int> &a, const pair<int,int> &b)-> bool {
            return a.second > b.second;
        });
        sorted_degree.clear();
        for(pair<int,int> p: tempMap)
            sorted_degree.push_back(p.first);
    }

    void setSContri(string path){
        s_contri = new double[vertexNum];
        s_contri_order = new int[vertexNum];
        memset(s_contri, 0, vertexNum);
        char line[256];
        int index = 0, node;
        double value;
        FILE *fd = fopen(path.c_str(), "r");
        vector<string> inStr;
        while(NULL != fgets(line, 256, fd)){
            boost::split(inStr, line, boost::is_any_of(" "));
            s_contri[node] = stod(inStr[1]);
            s_contri_order[index] = stoi(inStr[0]);
            index++;
        }
        fclose(fd);
        is_s_contri = true;
        s_contri_path = path;
    }

    void showSContri(){
        for(int i = 0;i < vertexNum;i++)
            printf("%d %lf\n", i, s_contri[i]);
    }

    void showData(){
        for(int i = 0;i < vertexNum;i++){
            for(int j = 0;j < outDegree[i];j++){
                printf("%d - %d\n", i, getNeighbor(i,j));
            }
        }
        printf("Edge number: %d\n", edgeNum);
    }

    double getProb(int cseed, int cseede) const {
        if(type == "IC")
            return IC_prob;
        else if(type == "VIC"){
            int index = [this, cseed, cseede]()-> int {
                for(int i = 0;i < outDegree[cseed];i++)
                    if(getNeighbor(cseed, i) == cseede)
                        return i;
                return -1;
            }();
            if(index < 0){
                cout << "error occurs in getProb" << endl;
                exit(1);
            }
            return getProbability(cseed, index);
        }
        else if(type == "WC")
            return 1./inDegree[cseede];
        else if(type == "LT")
            return 0.;
        else
            printf("Invalid model\n");
            return 0.;
    }
};
