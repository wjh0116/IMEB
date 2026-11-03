#ifndef __GREEDY_H__
#define __GREEDY_H__

#include "head.h"
#include "iheap.h"

struct GreedyResult {
    vector<Edge> B;
    double Lambda = 0.0;
    double Lambda_u = 0.0;
};

class Pool {
public:
    long long n_nodes;
    unordered_map<long long, int> eid;
    vector<Edge> edges;
    vector<vector<int>> edge_samples;
    vector<vector<int>> sample_edges;
    vector<double> weight;
    int num_samples = 0;

    explicit Pool(long long n) : n_nodes(n) {}

    inline long long key(int u, int w) const { return (long long)u * n_nodes + w; }

    int getEid(int u, int w) {
        long long kk = key(u, w);
        auto it = eid.find(kk);
        if (it != eid.end()) return it->second;
        int id = (int)edges.size();
        eid[kk] = id;
        edges.push_back(Edge(u, w));
        edge_samples.emplace_back();
        return id;
    }

    void addSample(const vector<Edge>& Q, double N) {
        int sid = num_samples++;
        weight.push_back(N);
        vector<int> es;
        es.reserve(Q.size());
        for (auto& e : Q) es.push_back(getEid(e.first, e.second));
        sort(es.begin(), es.end());
        es.erase(unique(es.begin(), es.end()), es.end());
        for (int id : es) edge_samples[id].push_back(sid);
        sample_edges.push_back(std::move(es));
    }

    void clear() {
        eid.clear(); edges.clear(); edge_samples.clear();
        sample_edges.clear(); weight.clear(); num_samples = 0;
    }
};

inline GreedyResult weightedGreedy(Pool& P, int k) {
    GreedyResult res;
    int E = (int)P.edges.size();
    if (E == 0) return res;

    vector<double> gain(E, 0.0);
    for (int e = 0; e < E; e++) {
        double g = 0.0;
        for (int s : P.edge_samples[e]) g += P.weight[s];
        gain[e] = g;
    }
    IHeap h; h.init(E);
    for (int e = 0; e < E; e++) h.insert(e, gain[e]);

    vector<char> covered(P.num_samples, 0);
    vector<char> removed(E, 0);
    double coverage = 0.0, Lu = 1e300;

    vector<double> tmp; tmp.reserve(E);
    auto topkSum = [&]() -> double {
        tmp.clear();
        for (int e = 0; e < E; e++) if (!removed[e] && gain[e] > 0) tmp.push_back(gain[e]);
        if ((int)tmp.size() > k) {
            nth_element(tmp.begin(), tmp.begin() + k, tmp.end(), greater<double>());
            tmp.resize(k);
        }
        double s = 0.0; for (double x : tmp) s += x; return s;
    };

    res.B.reserve(k);
    for (int i = 0; i <= k; i++) {

        Lu = min(Lu, coverage + topkSum());
        if (i == k || h.empty()) break;

        auto [e, g] = h.extractMax();
        if (g <= 0.0) break;
        removed[e] = 1;
        res.B.push_back(P.edges[e]);
        coverage += gain[e];

        for (int s : P.edge_samples[e]) {
            if (covered[s]) continue;
            covered[s] = 1;
            double w = P.weight[s];
            for (int e2 : P.sample_edges[s]) {
                if (removed[e2]) continue;
                gain[e2] -= w;
                h.updateKey(e2, gain[e2]);
            }
        }
    }
    res.Lambda = coverage;
    res.Lambda_u = Lu;
    return res;
}

inline double evalCoverage(Pool& P, const vector<Edge>& B) {
    vector<char> covered(P.num_samples, 0);
    double cov = 0.0;
    for (auto& e : B) {
        auto it = P.eid.find(P.key(e.first, e.second));
        if (it == P.eid.end()) continue;
        for (int s : P.edge_samples[it->second]) {
            if (!covered[s]) { covered[s] = 1; cov += P.weight[s]; }
        }
    }
    return cov;
}

#endif
