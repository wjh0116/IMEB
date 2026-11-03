#ifndef __CRITICALEDGE_H__
#define __CRITICALEDGE_H__

#include "graph.h"
#include "DESCE.h"

class CriticalEdge {
public:
    Graph& G;
    vector<int> seeds;
    int theta;

    map<Edge, long long> count;

    vector<bool> is_seed;
    SFMTRng gen;

    CriticalEdge(Graph& G, vector<int>& seeds, int theta)
        : G(G), seeds(seeds), theta(theta), is_seed(G.n, false) {
        unsigned s = chrono::system_clock::now().time_since_epoch().count();
        gen.seed((uint32_t)s);
        for (int s : seeds) is_seed[s] = true;
    }

    void run() {
        auto t_start = chrono::high_resolution_clock::now();
        count.clear();

        for (int t = 0; t < theta; t++) {
            if ((t + 1) % 1000 == 0 || t == 0) {
                cout << "CriticalEdge: sample " << (t + 1) << "/" << theta << endl;
            }
            processSample();
        }

        auto t_end = chrono::high_resolution_clock::now();
        double elapsed = chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count() / 1000.0;
        cout << "CriticalEdge completed in " << elapsed << " seconds (" << theta << " samples)" << endl;
    }

private:
    void processSample() {

        vector<bool> forward_visited(G.n, false);
        queue<int> fwd_q;

        vector<vector<int>> sampled_adj(G.n);
        vector<vector<int>> sampled_rev(G.n);
        vector<int> reachable;

        for (int s : seeds) {
            if (!forward_visited[s]) {
                forward_visited[s] = true;
                fwd_q.push(s);
            }
        }

        while (!fwd_q.empty()) {
            int u = fwd_q.front();
            fwd_q.pop();
            for (size_t j = 0; j < G.adj[u].size(); j++) {
                if (G.blocked[u][j]) continue;
                int w = G.adj[u][j];
                double p = G.prob[u][j];
                if (gen.randomDouble() < p) {
                    sampled_adj[u].push_back(w);
                    sampled_rev[w].push_back(u);
                    if (!forward_visited[w]) {
                        forward_visited[w] = true;
                        fwd_q.push(w);
                        if (!is_seed[w]) {
                            reachable.push_back(w);
                        }
                    }
                }
            }
        }

        if (reachable.empty()) return;

        int v = reachable[gen.randomUint32() % reachable.size()];

        vector<int> node_map(G.n, -1);
        vector<int> h_nodes;

        node_map[v] = h_nodes.size();
        h_nodes.push_back(v);

        queue<int> bwd_q;
        bwd_q.push(v);
        while (!bwd_q.empty()) {
            int node = bwd_q.front();
            bwd_q.pop();
            for (int u : sampled_rev[node]) {
                if (node_map[u] != -1) continue;
                if (!forward_visited[u]) continue;
                node_map[u] = h_nodes.size();
                h_nodes.push_back(u);
                bwd_q.push(u);
            }
        }

        int nh = h_nodes.size();
        if (nh <= 1) return;

        struct HEdge { int u_orig, w_orig; int u_h, w_h; };
        vector<HEdge> h_edges;

        for (int ui = 0; ui < nh; ui++) {
            int u = h_nodes[ui];
            for (int w : sampled_adj[u]) {
                if (node_map[w] == -1) continue;
                int wi = node_map[w];
                h_edges.push_back({u, w, ui, wi});
            }
        }

        int super_src = nh;
        DomTree dt;
        dt.init(nh + 1, super_src);

        for (int s : seeds) {
            if (node_map[s] != -1) {
                dt.addEdge(super_src, node_map[s]);
            }
        }

        for (auto& e : h_edges) {
            dt.addEdge(e.u_h, e.w_h);
        }
        vector<int>& idom = dt.computeIdom();

        int v_h = node_map[v];
        vector<bool> dominates_v(nh + 1, false);
        {
            int cur = v_h;
            while (cur != -1 && !dominates_v[cur]) {
                dominates_v[cur] = true;
                if (cur == idom[cur]) break;
                cur = idom[cur];
            }
        }

        vector<vector<int>> dom_children(nh + 1);
        for (int i = 0; i < nh; i++) {
            if (idom[i] != -1 && idom[i] != i) {
                dom_children[idom[i]].push_back(i);
            }
        }

        vector<int> dfs_in(nh + 1, -1), dfs_out(nh + 1, -1);
        {
            int timer = 0;
            vector<pair<int,int>> stk;
            stk.push_back({super_src, 0});
            dfs_in[super_src] = timer++;
            while (!stk.empty()) {
                auto& [nd, ci] = stk.back();
                if (ci < (int)dom_children[nd].size()) {
                    int c = dom_children[nd][ci++];
                    dfs_in[c] = timer++;
                    stk.push_back({c, 0});
                } else {
                    dfs_out[nd] = timer++;
                    stk.pop_back();
                }
            }
        }

        auto dom_dominates = [&](int a, int b) -> bool {
            if (dfs_in[a] == -1 || dfs_in[b] == -1) return false;
            return dfs_in[a] <= dfs_in[b] && dfs_out[b] <= dfs_out[a];
        };

        vector<vector<int>> h_in_adj(nh);
        for (int s : seeds) {
            if (node_map[s] != -1) {
                h_in_adj[node_map[s]].push_back(super_src);
            }
        }
        for (auto& e : h_edges) {
            h_in_adj[e.w_h].push_back(e.u_h);
        }

        long long Nw = (long long)reachable.size();
        for (int wi = 0; wi < nh; wi++) {
            if (!dominates_v[wi]) continue;
            int cw = 0, ustar = -1;
            for (int wp : h_in_adj[wi]) {
                if (!dom_dominates(wi, wp)) { ustar = wp; if (++cw > 1) break; }
            }
            if (cw == 1 && ustar < nh)
                count[{h_nodes[ustar], h_nodes[wi]}] += Nw;
        }
    }
};

#endif
