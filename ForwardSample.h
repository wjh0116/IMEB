#ifndef __FORWARD_SAMPLE_H__
#define __FORWARD_SAMPLE_H__

#include "graph.h"

class ForwardSampler {
public:
    Graph& G;
    vector<int> seeds;
    vector<char> is_seed;
    SFMTRng& gen;

    vector<char>        visited;
    vector<int>         dist;
    vector<vector<int>> sadj;
    vector<vector<int>> srev;
    vector<int>         reachable;
    vector<int>         touched;
    vector<int>         bfsq;
    int v = -1;
    int N = 0;

    ForwardSampler(Graph& G_, const vector<int>& seeds_, SFMTRng& gen_)
        : G(G_), seeds(seeds_), is_seed(G_.n, 0), gen(gen_),
          visited(G_.n, 0), dist(G_.n, -1), sadj(G_.n), srev(G_.n) {
        for (int s : seeds) is_seed[s] = 1;
    }

    bool sample() {

        for (int x : touched) { visited[x] = 0; dist[x] = -1; sadj[x].clear(); srev[x].clear(); }
        touched.clear();
        reachable.clear();

        bfsq.clear();
        for (int s : seeds) {
            if (!visited[s]) {
                visited[s] = 1; dist[s] = 0;
                touched.push_back(s);
                bfsq.push_back(s);
            }
        }

        size_t head = 0;
        while (head < bfsq.size()) {
            int u = bfsq[head++];
            const auto& au = G.adj[u];
            const auto& pu = G.prob[u];
            const auto& bu = G.blocked[u];
            for (size_t j = 0; j < au.size(); j++) {
                if (bu[j]) continue;
                if (gen.randomDouble() < pu[j]) {
                    int w = au[j];
                    sadj[u].push_back(w);
                    srev[w].push_back(u);
                    if (!visited[w]) {
                        visited[w] = 1; dist[w] = dist[u] + 1;
                        touched.push_back(w);
                        bfsq.push_back(w);
                        if (!is_seed[w]) reachable.push_back(w);
                    }
                }
            }
        }

        if (reachable.empty()) { N = 0; v = -1; return false; }
        v = reachable[gen.randomUint32() % reachable.size()];
        N = (int)reachable.size();
        return true;
    }
};

#endif
