#pragma once

#include <cstdio>
#include <cstdlib>
#include <set>
#include <map>
#include <vector>
#include <bitset>
#include <string>
#include <random>
#include <iostream>

#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>
#include <boost/thread.hpp>
#include <boost/bind/bind.hpp>
#include <boost/thread/mutex.hpp>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "network.hpp"
#include "rtuple.hpp"
#include "tools.hpp"
#include "macro.hpp"

class DiffusionState_MIC{
    short *temp_state_1, *temp_state_2;
    std::string pri_type, type_bak;
    std::vector<std::mt19937> randn;
    std::set<int> new_active[THREAD], new_active_temp[THREAD];
    boost::mutex mt;

    void readPriority(int num){
        std::string path = "../data/perms/"+std::to_string(num)+"_perms.txt";
        char line[50];
        int index = 0;
        FILE *fd = fopen(path.c_str(), "r");
        std::vector<std::vector<int>> temp;
        while((NULL != fgets(line, 50, fd)) && index < _max){
            std::vector<std::string> inStr;
            std::vector<int> tempp;
            boost::split(inStr, line, boost::is_any_of(" "));
            for(std::string str: inStr){
                if(str == "\r\n" || str == "\n")    break;
                tempp.push_back(std::stoi(str));
            }
            temp.push_back(tempp);
            index++;
        }
        caspriority[num] = temp;
        fclose(fd);
    }

    int pri_rand(const std::set<std::pair<int,short>> &pri, int tid){
        int n = randn[tid]() % pri.size();
        for(std::pair<int,short> p: pri){
            switch(pri_type[0]){
            case 'u':
                if(p.second == 4)
                    return 4;
                break;
            case 'l':
                if(p.second == 4)
                    n++;
                break;
            }
            if(n == 0)
                return p.second;
            --n;
        }
        return pri.begin()->second;
    }

    int pri_cas(const std::set<std::pair<int,short>> &pri, int node, int shift){
        std::vector<int> prior = caspriority[10][node%_max];
        int maxpri = pri.begin()->second;
        for(std::pair<int,short> p: pri){
            switch(pri_type[0]){
            case 'u':
                if(p.second == 4)
                    return 4;
            case 'l':
                if(p.second == 4)
                    continue;
            }
            if(prior[p.second] > prior[maxpri])
                maxpri = p.second;
        }
        return maxpri;
    }

    int pri_cas(const std::set<short> &pri, int node, int shift){
        std::vector<int> prior = caspriority[10][node%_max];
        int maxpri = *(pri.begin());
        for(short p: pri){
            switch(pri_type[0]){
            case 'u':
                if(p == 4)
                    return 4;
            case 'l':
                if(p == 4)
                    continue;
            }
            if(prior[p] > prior[maxpri])
                maxpri = p;
        }
        return maxpri;
    }

    int pri_nei(const std::set<std::pair<int,short>> &pri){
        int maxpri = -1, minnode = vnum;
        for(std::pair<int,short> p: pri){
            switch(pri_type[0]){
            case 'u':
                if(p.second == 4)
                    return 4;
            case 'l':
                if(p.second == 4)
                    continue;
            }
            if(p.first < minnode){
                minnode = p.first;
                maxpri = p.second;
            }
        }
        if(maxpri == -1)
            return pri.begin()->second;
        return maxpri;
    }

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

    int priority(const std::set<short> &pri, int n, int shift=0){
        return pri_cas(pri, n, shift);
    }

    void diffuseOneRound(const Network &network, short *state, int tid){
        int cseede, base = vnum * tid;
        double prob, rd;
        new_active_temp[tid].clear();
        std::map<int,std::set<std::pair<int,short>>> cas_pri;
        for(int cseed: new_active[tid]){
            for(int i = 0;i < network.outDegree[cseed];i++){
                cseede = network.getNeighbor(cseed, i);
                prob = network.getProb(cseed, cseede);
                rd = (double)randn[tid]()/randn[tid].max();
                if(rd < prob && state[cseede] == -1){
                    if(cas_pri.find(cseede) == cas_pri.end())
                        cas_pri[cseede] = std::set<std::pair<int,short>>();
                    cas_pri[cseede].insert(std::make_pair(cseed, state[cseed]));
                    if(new_active_temp[tid].find(cseede) == new_active_temp[tid].end())
                        new_active_temp[tid].insert(cseede);
                }
            }
        }
        new_active[tid].clear();
        for(int i: new_active_temp[tid]){
            new_active[tid].insert(i);
            state[i] = priority(cas_pri[i], i, tid);
        }
    }

    void diffuseOneRound(short *state, rTuple &rtup, int tid){
        int base = tid * vnum;
        new_active_temp[tid].clear();
        std::map<int,std::set<std::pair<int,short>>> cas_pri;
        for(int cseed: new_active[tid]){
            for(int cseede: rtup.relations[cseed]){
                if(state[cseede] == -1){
                    if(cas_pri.find(cseede) == cas_pri.end())
                        cas_pri[cseede] = std::set<std::pair<int,short>>();
                    cas_pri[cseede].insert(std::make_pair(cseed, state[cseed]));
                    if(new_active_temp[tid].find(cseede) == new_active_temp[tid].end())
                        new_active_temp[tid].insert(cseede);
                }
            }
        }
        new_active[tid].clear();
        for(int i: new_active_temp[tid]){
            new_active[tid].insert(i);
            state[i] = priority(cas_pri[i], i, tid, 1);
        }
    }

    int reSpreadOneRound(const Network &network, short *state, rTuple &rtup, int tid){
        int cseede;
        double prob, rd;
        new_active_temp[tid].clear();
        for(int cseed: new_active[tid]){
            for(int j = 0;j < network.inDegree[cseed];j++){
                cseede = network.getReverseNeighbor(cseed, j);
                prob = network.getProb(cseede, cseed);
                if((double)randn[tid]()/randn[tid].max() < prob && state[cseede] != -2){
                    if(rtup.relations.find(cseede) == rtup.relations.end()){
                        rtup.relations[cseede] = std::set<int>();
                        new_active_temp[tid].insert(cseede);
                    }
                    rtup.relations[cseede].insert(cseed);
                }
            }
        }
        new_active[tid].clear();
        bool islast = false;
        for(int i: new_active_temp[tid]){
            if(state[i] >= 0){
                islast = true;
                rtup.seed.insert(i);
            } else {
                new_active[tid].insert(i);
                state[i] = -2;
            }
            rtup.upper.insert(i);
        }
        if(islast){
            new_active[tid].clear();
            rtup.isdiff = true;
        } else
            for(int i: new_active_temp[tid])
                rtup.lower.insert(i);
        return 0;
    }

    int reSpreadOnce(const Network &network, int cindex, rTuple &rtup, int tid){
        int base = vnum * tid;
        for(int i = 0;i < vnum;i++) temp_state_1[base+i] = seed_state[i];
        temp_state_1[base+cindex] = -2;
        new_active[tid].clear();
        rtup.lower.insert(cindex);
        new_active[tid].insert(cindex);
        while(!new_active[tid].empty())
            reSpreadOneRound(network, temp_state_1+base, rtup, tid);
        return 0;
    }

    double getRTuple(const Network &network, rTuple &rtup){
        int tid = -1;
        mt.lock();
        tid = scheduler._Find_first();
        if(tid >= THREAD)
            std::cout << "get r tuple" << std::endl;
        scheduler.flip(tid);
        mt.unlock();
        int cindex = randn[tid]()%vnum;
        rtup.clear();
        rtup.node_v = cindex;
        rtup.upper.insert(cindex);
        rtup.relations[cindex] = std::set<int>();
        if(seed_state[cindex] != -1){
            rtup.seed.insert(cindex);
            mt.lock();
            scheduler.flip(tid);
            mt.unlock();
            return 0;
        }
        int ret = reSpreadOnce(network, cindex, rtup, tid);
        mt.lock();
        scheduler.flip(tid);
        mt.unlock();
        return ret;
    }

    void updateCascade(){
        memset(seed_state, -1, vnum*sizeof(short));
        std::map<int,std::set<short>> seed_cascade;     // node, cascade
        for(std::pair<int,std::set<int>> p: seedsets){  // cascade, seeds
            for(int n: p.second){
                if(seed_cascade.find(n) == seed_cascade.end())
                    seed_cascade[n] = std::set<short>();
                seed_cascade[n].insert(p.first);
            }
        }
        for(std::pair<int,std::set<short>> p: seed_cascade)
            seed_state[p.first] = priority(p.second, p.first);
    }

    template <typename T>
    void freeSpace(T *p){
        if(p)   delete [] p;
    }

public:
    short *seed_state;
    int cnum, vnum, _max;
    std::set<int> seednodes;
    std::map<int,std::set<int>> seedsets;
    std::map<int,std::vector<std::vector<int>>> caspriority, caspriority_upper, caspriority_lower;
    std::bitset<THREAD> scheduler;

    DiffusionState_MIC(const Network &network, const std::string pri_type, mt19937 &mrand){
        cnum = 0;
        vnum = network.vertexNum;
        _max = 13;
        seed_state = nullptr;
        temp_state_1 = nullptr;
        temp_state_2 = nullptr;
        seed_state = new short[vnum];
        temp_state_1 = new short[vnum*THREAD];
        temp_state_2 = new short[vnum*THREAD];
        this->pri_type = pri_type;
        type_bak = pri_type;
        for(int i = 0;i < THREAD;i++)
            randn.push_back(mt19937(mrand()));
        memset(seed_state, -1, vnum*sizeof(short));
        
        readPriority(1);
        readPriority(3);
        readPriority(5);
        readPriority(10);
    }

    DiffusionState_MIC(const DiffusionState_MIC &diffusionState){
        cnum = diffusionState.cnum;
        _max = diffusionState._max;
        vnum = diffusionState.vnum;
        seed_state = nullptr;
        temp_state_1 = nullptr;
        temp_state_2 = nullptr;
        seed_state = new short[vnum];
        memcpy(seed_state, diffusionState.seed_state, vnum*sizeof(short));
        temp_state_1 = new short[vnum*THREAD];
        temp_state_2 = new short[vnum*THREAD];
        seednodes = diffusionState.seednodes;
        seedsets = diffusionState.seedsets;
        caspriority = diffusionState.caspriority;
        pri_type = diffusionState.pri_type;
        type_bak = pri_type;
        for(int i = 0;i < THREAD;i++)
            randn.push_back(mt19937(diffusionState.randn[i]));
    }

    ~DiffusionState_MIC(){
        freeSpace(seed_state);
        freeSpace(temp_state_1);
        freeSpace(temp_state_2);
    }

    void diffuse(const Network &network, int *result, int cindex, int round, int j){
        // auto start = std::chrono::high_resolution_clock::now();
        int tid = -1;
        mt.lock();
        tid = scheduler._Find_first();
        scheduler.flip(tid);
        mt.unlock();
        new_active[tid] = seednodes;
        int base = vnum * tid;
        for(int i = 0;i < vnum;i++) temp_state_1[base+i] = seed_state[i];
        for(int i = 0;i < round;i++){
            diffuseOneRound(network, temp_state_1+base, tid);
            if(new_active[tid].empty())  break;
        }
        for(int i = 0;i < vnum;i++)
            if(temp_state_1[base+i] == cindex)
                result[j]++;
        mt.lock();
        scheduler.flip(tid);
        mt.unlock();
        // auto end = std::chrono::high_resolution_clock::now();
        // printTime(start, end);
    }

    int seed(const std::set<int> &seed_set){
        std::set<int> new_seed;
        for(int i: seed_set){
            seednodes.insert(i);
            new_seed.insert(i);
        }
        seedsets[cnum++] = new_seed;
        updateCascade();
        return cnum-1;
    }

    void removeSeed(int cindex){
        if(seedsets.find(cindex) == seedsets.end()){
            std::cout << "WARNING: removing cascade does not exist" << std::endl;
            return;
        }
        for(int i: seedsets[cindex]){
            // seed_state[i] = -1;
            seednodes.erase(i);
        }
        seedsets.erase(cindex);
        cnum--;
        updateCascade();
    }

    double getRTuples(const Network &network, std::vector<rTuple> &rtup, double size){
        auto start = std::chrono::high_resolution_clock::now();
        int countdiff = 0;
        int rtup_size = rtup.size();
        if(rtup.empty())    rtup = std::vector<rTuple>(size);
        else    rtup.resize(size+rtup.size());
        int new_size = rtup.size();
        boost::asio::thread_pool pool(THREAD);
        allSet(scheduler);
        for(int i = rtup_size;i < new_size;i++){
            auto bind_fn = boost::bind(&DiffusionState_MIC::getRTuple, this, ref(network), ref(rtup[i]));
            boost::asio::post(pool, bind_fn);
        }
        pool.join();
        for(int i = rtup_size;i < new_size;i++)
            if(rtup[i].isdiff)
                countdiff++;
        auto end = std::chrono::high_resolution_clock::now();
        // printTime(start,end, true,"get r tuple");
        return (double)countdiff;
    }

    double expInfluenceComplete_new(const Network &network, std::vector<rTuple> &rtup, std::set<int> &solution){
        return computeG(solution, rtup, vnum, "mid", nullptr);
    }

    double expInfluenceComplete(const Network &network, int times, int cindex){
        int c_result[times], tid = 0;
        double result = 0.;
        memset(c_result, 0, times*sizeof(int));
        boost::asio::thread_pool pool(THREAD);
        allSet(scheduler);
        for(int i = 0;i < times;i++){
            auto bind_fn = boost::bind(&DiffusionState_MIC::diffuse, this, ref(network), c_result, cindex, vnum, tid);
            boost::asio::post(pool, bind_fn);
            tid++;
        }
        pool.join();
        for(int j = 0;j < times;j++)
            result += c_result[j];
        return result/times;
    }

    bool computeMid_g(const std::set<int> &seed, rTuple &rtup, int tid){
        new_active[tid].clear();
        int base = vnum * tid, cid;
        for(int i: rtup.upper)
            temp_state_1[base+i] = -1;
        for(int i: rtup.seed){
            temp_state_1[base+i] = seed_state[i];
            if(seed_state[i] > -1)
                new_active[tid].insert(i);
            else{
                std::cout << "ERROR: seed state" << std::endl;
                exit(1);
            }
        }
        for(int i: seed){
            if(rtup.upper.find(i) != rtup.upper.end()){
                cid = cnum;
                if(seed_state[i] > -1){
                    std::set<std::pair<int,short>> temp_set;
                    temp_set.insert(std::make_pair(0, seed_state[i]));
                    temp_set.insert(std::make_pair(0, cnum));
                    cid = priority(temp_set, i, tid, 0);
                }
                temp_state_1[base+i] = cid;
                new_active[tid].insert(i);
            }
        }
        for(int i = 0;i < 2*rtup.relations.size();i++){
            if(new_active[tid].find(rtup.node_v) != new_active[tid].end()){
                if(temp_state_1[base+rtup.node_v] == cnum)   return true;
                return false;
            } else if(new_active[tid].empty())
                return false;
            diffuseOneRound(temp_state_1+base, rtup, tid);
        }
        std::cout << "ERROR: unusual exit from DiffusionState_MIC.compute_g" << std::endl;
        exit(1);
    }

    int compute_g(const std::set<int> &seed, rTuple &rtup, const std::string &type, int *result, int tid){
        int temp = 0;
        bool flag = false;
        if(tid < 0){
            flag = true;
            mt.lock();
            tid = scheduler._Find_first();
            scheduler.flip(tid);
            mt.unlock();
        }
        switch(type[0]){
        case 'u':
            if(intersection(seed, rtup.upper))  temp = 1;
            else    temp = -1;
            break;
        case 'm':
            if(intersection(seed, rtup.lower))  temp = 2;
            else if(intersection(seed, rtup.upper)){
                if(computeMid_g(seed, rtup, tid))    temp = 1;
                else    temp = -1;
            } else  temp = -2;
            break;
        case 'l':
            if(intersection(seed, rtup.lower))  temp = 2;
            else    temp = -2;
            break;
        }
        if(flag){
            mt.lock();
            scheduler.flip(tid);
            mt.unlock();
        }
        if(result)  *result = temp;
        return temp;
    }

    void __parallel(std::set<int> &S, std::vector<rTuple>::iterator _first, std::vector<rTuple>::iterator _second, const std::string type, int *result, int tid){
        auto start = std::chrono::high_resolution_clock::now();
        int count = 0;
        for(auto it = _first;it != _second;it++)
            if(compute_g(S, *it, type, nullptr, tid) > 0)
                count++;
        *result = count;
        auto end = std::chrono::high_resolution_clock::now();
        // printTime(start, end, true, std::to_string(tid));
    }

    double computeG(std::set<int> &S, std::vector<rTuple> &rtup, int n, const std::string &type, double *result){
        int count = 0, each_thread = rtup.size()/THREAD, rest = rtup.size()-THREAD*each_thread;
        double output;
        int results[THREAD];
        memset(results, 0, sizeof(results));
        std::vector<boost::thread> thread_list;
        std::vector<rTuple>::iterator head = rtup.begin(), tail;
        auto start = std::chrono::high_resolution_clock::now();
        for(int i = 0;i < THREAD;i++){
            tail = head+each_thread;
            if(rest){
                tail++;
                rest--;
            }
            if(tail > rtup.end())   tail = rtup.end();
            boost::thread new_thread(
                boost::bind(&DiffusionState_MIC::__parallel, this, ref(S), head, tail, ref(type), results+i, i)
            );
            head = tail;
            thread_list.push_back(boost::move(new_thread));
        }
        for(boost::thread &t: thread_list)
            t.join();
        auto end = std::chrono::high_resolution_clock::now();
        // printTime(start, end, true, "total");
        // for(rTuple &rt: rtup){
        //     auto bind_fn = boost::bind(&DiffusionState_MIC::compute_g, this, ref(S), ref(rt), ref(type), results+tid, -1);
        //     boost::asio::post(pool, bind_fn);
        //     tid++;
        // }
        // pool.join();
        for(int i = 0;i < THREAD;i++)
            count += results[i];
        output = (double)n*count/rtup.size();
        if(result)  *result = output;
        return output;
    }

    double computeG_old(std::set<int> &S, std::vector<rTuple> &rtup, int n, const std::string &type, double *result){
        int count = 0, tid = 0;
        double output;
        int *results = new int[rtup.size()];
        boost::asio::thread_pool pool(THREAD);
        allSet(scheduler);
        for(rTuple &rt: rtup){
            auto bind_fn = boost::bind(&DiffusionState_MIC::compute_g, this, ref(S), ref(rt), ref(type), results+tid, -1);
            boost::asio::post(pool, bind_fn);
            tid++;
        }
        pool.join();
        for(int i = 0;i < rtup.size();i++)
            if(results[i] > 0)
                count++;
        output = (double)n*count/rtup.size();
        if(result)  *result = output;
        return output;
    }

    void setBack(){
        pri_type = type_bak;
    }

    void setUpper(){
        type_bak = pri_type;
        pri_type = "upper";
    }

    void setLower(){
        type_bak = pri_type;
        pri_type = "lower";
    }
};
