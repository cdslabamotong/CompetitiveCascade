#pragma once

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;

class sortedMap{
    int _size;
public:
    vector<int> node_order;
    map<int,double> value_list;
    sortedMap(){
        _size = 0;
    }
    sortedMap(const sortedMap& sm){
        node_order = sm.node_order;
        value_list = sm.value_list;
        _size = sm._size;
    }
    inline int size() const {
        return _size;
    }
    void push_back(int n, double v){
        node_order.push_back(n);
        value_list[n] = v;
        _size++;
    }
    void remove(int n){
        vector<int>::iterator v_iter;
        if((v_iter = find(node_order.begin(), node_order.end(), n)) == node_order.end())
            cout << "remove fail: no node exists" << endl;
        else{
            node_order.erase(v_iter);
            if(value_list.find(n) == value_list.end())
                cout << "remove fail: node is associated with null" << endl;
            else{
                value_list.erase(n);
                _size--;
            }
        }
    }
    int get(int id){
        if(id > node_order.size())
            cout << "sortedMap get out of bound" << endl;
        return node_order[id];
    }
    double getvalue(int n){
        if(value_list.find(n) == value_list.end())
            cout << "sortedMap getvalue node not exist" << endl;
        return value_list[n];
    }
    bool contains(int n){
        return find(node_order.begin(), node_order.end(), n) != node_order.end();
    }
    int insert(int n, double v, int l, int r){
        int mid = -1;
        if(v >= value_list[node_order[l]]){
            node_order.insert(node_order.begin()+l, n);
            _size++;
            return l;
        } else if (v <= value_list[node_order[r-1]]){
            node_order.insert(node_order.begin()+r, n);
            _size++;
            return r;
        } else if(r-l < 10){
            for(int i = l;i < r;i++){
                if(v >= value_list[node_order[i]]){
                    node_order.insert(node_order.begin()+i, n);
                    _size++;
                    return i;
                }
            }
            cout << "something wrong when insert" << endl;
        } else{
            mid = (l+r)>>1;
            if(v >= value_list[node_order[mid]])
                return insert(n, v, l, mid);
            return insert(n, v, mid, r);
        }
    }
    int update(int n, double v){
        value_list[n] = v;
        vector<int>::iterator v_iter = find(node_order.begin(), node_order.end(), n);
        if(v_iter == node_order.end()){
            cout << n << " does not in the list" << endl;
            return -1;
        }
        else{
            node_order.erase(v_iter);
            return insert(n, v, 0, node_order.size());
        }
    }
};
