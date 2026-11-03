#ifndef __INFGRAPH_H__
#define __INFGRAPH_H__

#include "graph.h"

class InfGraph : public Graph {
public:

    vector<int> seed_set;

    SFMTRng gen;

    InfGraph(string folder, string graph_file)
        : Graph(folder, graph_file) {
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        gen.seed((uint32_t)seed);
    }

    void setSeedSet(const vector<int>& seeds) {
        seed_set = seeds;
    }

    double estimateInfluence(int num_simulations) {
        double total = 0.0;
        for (int sim = 0; sim < num_simulations; sim++) {
            total += simulateIC();
        }
        return total / num_simulations;
    }

    int simulateIC() {
        vector<bool> activated(n, false);
        queue<int> q;
        for (int s : seed_set) {
            activated[s] = true;
            q.push(s);
        }

        int count = seed_set.size();
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (size_t i = 0; i < adj[u].size(); i++) {
                int v = adj[u][i];
                if (activated[v]) continue;

                double p = getEffectiveProb(u, i);
                if (gen.randomDouble() < p) {
                    activated[v] = true;
                    q.push(v);
                    count++;
                }
            }
        }
        return count;
    }

};

#endif
