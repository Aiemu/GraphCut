#ifndef GRAPHCUT_MAXFLOW_HPP
#define GRAPHCUT_MAXFLOW_HPP

#define EPS 1e-8
#define INF 1e10

#include <string>
#include <cstring>
#include <queue>
#include <cmath>
#include <algorithm>

// Max Flow Algorithm
class MaxFlow{
public:
    MaxFlow(int, int);
    ~MaxFlow();

    void search(bool*);
    void insert(int, int, double);
    void reset();
    double dinitz();

    int get_begin();
    int get_end();
    int get_n_num();

    void set_begin(int);
    void set_end(int);
    void set_n_num(int);

private:
    int begin, end, num, n_num;
    int *head, *next, *pos, *node_list;
    double *res;

    std::queue<int> que;

    void insert_d(int, int, double);
    bool searchBySide();
    double searchByPath(int, double);
};

MaxFlow::MaxFlow(int len_q, int len_r) {
    begin = 0;
    end = num = 1;
    n_num = len_q;

    head = new int[len_q];
    node_list = new int[len_q];

    next = new int[len_r];
    pos = new int[len_r];
    res = new double[len_r];

    reset();
}
MaxFlow::~MaxFlow() = default;

int MaxFlow::get_begin() {
    return begin;
}

int MaxFlow::get_end() {
    return end;
}

int MaxFlow::get_n_num() {
    return n_num;
}

void MaxFlow::set_begin(int _begin) {
    begin = _begin;
}

void MaxFlow::set_end(int _end) {
    end = _end;
}

void MaxFlow::set_n_num(int _num) {
    n_num = _num;
}

void MaxFlow::reset() {
    num = 1;

    for (int i = 0; i < n_num + 1; i++) {
        head[i] = 0;
    }
}

double MaxFlow::dinitz() {
    double ret = 0;

    for (; searchBySide(); ){
        ret += searchByPath(begin, INF);
    }

    return ret;
}

void MaxFlow::insert(int i, int _pos, double _res) {
    insert_d(i, _pos, _res);
    insert_d(_pos, i, _res);
}

void MaxFlow::insert_d(int i, int _pos, double _res) {
    num++;

    res[num] = _res;
    pos[num] = _pos;
    next[num] = head[i];

    head[i] = num;
}

void MaxFlow::search(bool *search_list) {
    que = std::queue<int>();
    que.push(begin);

    for(int i = 0; i <= n_num; i++) {
        search_list[i] = false;
        node_list[i] = 0;
    }
    node_list[begin] = 1;

    for (; !que.empty(); ){
        int he = que.front();

        que.pop();
        for (int i = head[he]; i; ){
            if (res[i] > 0 && !node_list[pos[i]]){
                que.push(pos[i]);
                node_list[pos[i]] = 1;
                search_list[pos[i]] = true;
            }
            i = next[i];
        }
    }
}

bool MaxFlow::searchBySide() {
    que = std::queue<int>();
    que.push(begin);

    for (int i = 0; i < n_num + 1; i++) {
        node_list[i] = -1;
    }

    node_list[begin] = 0;

    for (; !que.empty(); ){
        int he = que.front();

        que.pop();

        for (int i = head[he]; i; ){
            if (node_list[pos[i]] == -1 && res[i] > 0){
                node_list[pos[i]] = node_list[he] + 1;

                que.push(pos[i]);
            }
            i = next[i];
        }
    }

    return node_list[end] != -1;
}

double MaxFlow::searchByPath(int _end, double _ret) {
    double delta, ret = 0;

    if (_end == end) {
        return _ret;
    }
    for (int i = head[_end]; i; ){
        if (node_list[pos[i]] == node_list[_end] + 1 && res[i] > 0){
            delta = _ret - ret;
            delta = searchByPath(pos[i], std::min(delta, res[i]));

            ret += delta;
            res[i] -= delta;
            res[i ^ 1] += delta;

            if (std::fabs(_ret - ret) < EPS) {
                return _ret;
            }
        }
        i = next[i];
    }

    if (std::fabs(ret) < EPS) {
        node_list[_end] = -1;
    }

    return ret;
}

#endif //GRAPHCUT_MAXFLOW_HPP
