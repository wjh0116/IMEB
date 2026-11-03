#ifndef __CRITICAL_EDGE_SET_H__
#define __CRITICAL_EDGE_SET_H__

#include "ForwardSample.h"
#include "DESCE.h"

class CESetBuilder {
public:
    Graph& G;
    DomTree dt;

    vector<int> node_map;
    vector<int> h_nodes;
    vector<int> bwd_q;

    explicit CESetBuilder(Graph& G_) : G(G_), node_map(G_.n, -1) {}

    void build(const ForwardSampler& fs, vector<Edge>& R) {
        int v = fs.v;

        h_nodes.clear();
        node_map[v] = 0;
        h_nodes.push_back(v);
        bwd_q.clear();
        bwd_q.push_back(v);
        for (size_t head = 0; head < bwd_q.size(); head++) {
            int node = bwd_q[head];
            for (int u : fs.srev[node]) {
                if (node_map[u] != -1) continue;
                node_map[u] = (int)h_nodes.size();
                h_nodes.push_back(u);
                bwd_q.push_back(u);
            }
        }

        int nh = (int)h_nodes.size();
        if (nh <= 1) { for (int x : h_nodes) node_map[x] = -1; return; }

        struct HEdge { int u_orig, w_orig, u_h, w_h; };
        vector<HEdge> h_edges;
        for (int ui = 0; ui < nh; ui++) {
            int u = h_nodes[ui];
            for (int w : fs.sadj[u]) {
                int wi = node_map[w];
                if (wi == -1) continue;
                h_edges.push_back({u, w, ui, wi});
            }
        }

        int super_src = nh;
        dt.init(nh + 1, super_src);
        for (int s : fs.seeds) {
            int si = node_map[s];
            if (si != -1) dt.addEdge(super_src, si);
        }
        for (auto& e : h_edges) dt.addEdge(e.u_h, e.w_h);
        vector<int>& idom = dt.computeIdom();

        int v_h = node_map[v];
        vector<char> dominates_v(nh + 1, 0);
        for (int cur = v_h; cur != -1 && !dominates_v[cur]; ) {
            dominates_v[cur] = 1;
            if (cur == idom[cur]) break;
            cur = idom[cur];
        }

        vector<vector<int>> dom_children(nh + 1);
        for (int i = 0; i < nh; i++)
            if (idom[i] != -1 && idom[i] != i) dom_children[idom[i]].push_back(i);
        vector<int> tin(nh + 1, -1), tout(nh + 1, -1);
        {
            int timer = 0;
            vector<pair<int,int>> stk;
            stk.push_back({super_src, 0});
            tin[super_src] = timer++;
            while (!stk.empty()) {
                auto& [nd, ci] = stk.back();
                if (ci < (int)dom_children[nd].size()) {
                    int c = dom_children[nd][ci++];
                    tin[c] = timer++;
                    stk.push_back({c, 0});
                } else { tout[nd] = timer++; stk.pop_back(); }
            }
        }
        auto dominates = [&](int a, int b) -> bool {
            if (tin[a] == -1 || tin[b] == -1) return false;
            return tin[a] <= tin[b] && tout[b] <= tout[a];
        };

        vector<vector<int>> h_in(nh);
        for (int s : fs.seeds) { int si = node_map[s]; if (si != -1) h_in[si].push_back(super_src); }
        for (auto& e : h_edges) h_in[e.w_h].push_back(e.u_h);

        for (int wi = 0; wi < nh; wi++) {
            if (!dominates_v[wi]) continue;
            int cw = 0, ustar = -1;
            for (int wp : h_in[wi]) {
                if (!dominates(wi, wp)) { ustar = wp; if (++cw > 1) break; }
            }
            if (cw == 1 && ustar < nh)
                R.push_back(Edge(h_nodes[ustar], h_nodes[wi]));
        }

        for (int x : h_nodes) node_map[x] = -1;
    }
};

#endif
