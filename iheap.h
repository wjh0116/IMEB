#ifndef __IHEAP_H__
#define __IHEAP_H__

#include <vector>
#include <algorithm>
using namespace std;

class IHeap {
public:
    vector<double> key;
    vector<int> heap;
    vector<int> pos;
    int size;

    IHeap() : size(0) {}

    void init(int n) {
        key.assign(n, 0.0);
        heap.resize(n);
        pos.assign(n, -1);
        size = 0;
    }

    bool empty() const { return size == 0; }

    void insert(int id, double k) {
        key[id] = k;
        heap[size] = id;
        pos[id] = size;
        size++;
        siftUp(pos[id]);
    }

    pair<int, double> extractMax() {
        int id = heap[0];
        double k = key[id];
        size--;
        heap[0] = heap[size];
        pos[heap[0]] = 0;
        pos[id] = -1;
        if (size > 0) siftDown(0);
        return make_pair(id, k);
    }


    void updateKey(int id, double k) {
        double old_k = key[id];
        key[id] = k;
        if (k > old_k) siftUp(pos[id]);
        else siftDown(pos[id]);
    }

private:
    void siftUp(int i) {
        while (i > 0) {
            int parent = (i - 1) / 2;
            if (key[heap[i]] > key[heap[parent]]) {
                swap(heap[i], heap[parent]);
                pos[heap[i]] = i;
                pos[heap[parent]] = parent;
                i = parent;
            } else break;
        }
    }

    void siftDown(int i) {
        while (2 * i + 1 < size) {
            int child = 2 * i + 1;
            if (child + 1 < size && key[heap[child + 1]] > key[heap[child]])
                child++;
            if (key[heap[child]] > key[heap[i]]) {
                swap(heap[i], heap[child]);
                pos[heap[i]] = i;
                pos[heap[child]] = child;
                i = child;
            } else break;
        }
    }
};

#endif
