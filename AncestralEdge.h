#ifndef __ANCESTRAL_EDGE_H__
#define __ANCESTRAL_EDGE_H__

#include "ForwardSample.h"

inline void buildAESet(const ForwardSampler& fs, vector<Edge>& P) {
    int w = fs.v;
    while (!fs.is_seed[w]) {
        const vector<int>& in = fs.srev[w];
        if (in.empty()) break;

        int best = in[0];
        for (size_t i = 1; i < in.size(); i++) {
            int u = in[i];
            if (fs.dist[u] < fs.dist[best] ||
                (fs.dist[u] == fs.dist[best] && u > best)) {
                best = u;
            }
        }
        P.push_back(Edge(best, w));
        w = best;
    }
}

#endif
