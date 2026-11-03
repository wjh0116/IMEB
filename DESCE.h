#ifndef __DESCE_H__
#define __DESCE_H__

#include "graph.h"

class DomTree {
public:
    int N;
    int root;
    vector<vector<int>> succ;
    vector<vector<int>> pred;

    vector<int> dfnum;
    vector<int> vertex;
    vector<int> parent;
    int dfs_count;

    vector<int> semi;
    vector<int> idom;
    vector<int> best;
    vector<int> ancestor;
    vector<vector<int>> bucket;

    DomTree() : N(0), root(0), dfs_count(0) {}

    void init(int n, int r) {
        N = n;
        root = r;
        succ.assign(N, {});
        pred.assign(N, {});
    }

    void addEdge(int u, int v) {
        succ[u].push_back(v);
        pred[v].push_back(u);
    }

    vector<int>& computeIdom() {
        dfnum.assign(N, -1);
        vertex.resize(N);
        parent.assign(N, -1);
        semi.resize(N);
        idom.assign(N, -1);
        best.resize(N);
        ancestor.assign(N, -1);
        bucket.assign(N, {});
        dfs_count = 0;

        dfs(root);

        for (int i = 0; i < dfs_count; i++) {
            int w = vertex[i];
            semi[w] = i;
            best[w] = w;
        }

        for (int i = dfs_count - 1; i >= 1; i--) {
            int w = vertex[i];

            for (int v : pred[w]) {
                if (dfnum[v] == -1) continue;
                int u = eval(v);
                if (semi[u] < semi[w]) {
                    semi[w] = semi[u];
                }
            }

            int semi_node = vertex[semi[w]];
            bucket[semi_node].push_back(w);

            link(parent[w], w);

            int p = parent[w];
            for (int v : bucket[p]) {
                int u = eval(v);
                idom[v] = (semi[u] == semi[v]) ? p : u;
            }
            bucket[p].clear();
        }

        for (int i = 1; i < dfs_count; i++) {
            int w = vertex[i];
            if (idom[w] != vertex[semi[w]]) {
                idom[w] = idom[idom[w]];
            }
        }

        idom[root] = root;
        return idom;
    }

private:
    void dfs(int r) {

        struct Frame { int node; int par; };
        vector<Frame> stack;
        stack.push_back({r, -1});

        while (!stack.empty()) {
            auto [v, p] = stack.back();
            stack.pop_back();
            if (dfnum[v] != -1) continue;

            dfnum[v] = dfs_count;
            vertex[dfs_count] = v;
            parent[v] = p;
            dfs_count++;

            for (int i = (int)succ[v].size() - 1; i >= 0; i--) {
                int w = succ[v][i];
                if (dfnum[w] == -1) {
                    stack.push_back({w, v});
                }
            }
        }
    }

    void link(int p, int w) {
        ancestor[w] = p;
    }

    int eval(int v) {
        if (ancestor[v] == -1) return v;
        compress(v);
        return best[v];
    }

    void compress(int v) {

        vector<int> path;
        int u = v;
        while (ancestor[u] != -1 && ancestor[ancestor[u]] != -1) {
            path.push_back(u);
            u = ancestor[u];
        }

        for (int i = (int)path.size() - 1; i >= 0; i--) {
            int w = path[i];
            if (semi[best[ancestor[w]]] < semi[best[w]]) {
                best[w] = best[ancestor[w]];
            }
            ancestor[w] = ancestor[ancestor[w]];
        }
    }
};

class DESCE {
public:
    Graph& G;
    vector<int> seeds;
    int theta;

    map<Edge, double> delta;

    SFMTRng gen;

    DESCE(Graph& G, vector<int>& seeds, int theta)
        : G(G), seeds(seeds), theta(theta) {
        unsigned s = chrono::system_clock::now().time_since_epoch().count();
        gen.seed((uint32_t)s);
    }

    void run() {
        auto t_start = chrono::high_resolution_clock::now();

        delta.clear();

        for (int t = 0; t < theta; t++) {
            if ((t + 1) % 100 == 0 || t == 0) {
                cout << "DESCE: sample " << (t + 1) << "/" << theta << endl;
            }
            processSample();
        }

        for (auto& [e, val] : delta) {
            val /= theta;
        }

        auto t_end = chrono::high_resolution_clock::now();
        double elapsed = chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count() / 1000.0;
        cout << "DESCE completed in " << elapsed << " seconds (" << theta << " samples)" << endl;
    }

private:
    void processSample() {

        struct SampledEdge { int u, v; };
        vector<SampledEdge> sampled_edges;

        vector<bool> visited(G.n, false);
        queue<int> bfs_q;
        for (int s : seeds) {
            if (!visited[s]) {
                visited[s] = true;
                bfs_q.push(s);
            }
        }
        while (!bfs_q.empty()) {
            int u = bfs_q.front();
            bfs_q.pop();
            for (size_t j = 0; j < G.adj[u].size(); j++) {
                if (G.blocked[u][j]) continue;
                int v = G.adj[u][j];
                double p = G.prob[u][j];
                if (gen.randomDouble() < p) {
                    sampled_edges.push_back({u, v});
                    if (!visited[v]) {
                        visited[v] = true;
                        bfs_q.push(v);
                    }
                }
            }
        }

        int n_real = G.n;
        int super_src = n_real;
        int n_virtual = sampled_edges.size();
        int total_nodes = n_real + 1 + n_virtual;

        DomTree dt;
        dt.init(total_nodes, super_src);

        for (int s : seeds) {
            dt.addEdge(super_src, s);
        }

        vector<Edge> virtual_to_edge(n_virtual);

        for (int i = 0; i < n_virtual; i++) {
            int u = sampled_edges[i].u;
            int v = sampled_edges[i].v;
            int w = n_real + 1 + i;

            dt.addEdge(u, w);
            dt.addEdge(w, v);

            virtual_to_edge[i] = {u, v};
        }

        vector<int>& idom = dt.computeIdom();

        vector<vector<int>> dom_children(total_nodes);
        for (int i = 0; i < total_nodes; i++) {
            if (idom[i] != -1 && idom[i] != i) {
                dom_children[idom[i]].push_back(i);
            }
        }

        vector<int> subtree_vv(total_nodes, 0);

        vector<int> order;
        order.reserve(total_nodes);
        {
            vector<bool> visited(total_nodes, false);
            struct Frame { int node; };
            vector<Frame> stk;
            stk.push_back({super_src});

            queue<int> bfs;
            bfs.push(super_src);
            visited[super_src] = true;
            while (!bfs.empty()) {
                int v = bfs.front(); bfs.pop();
                order.push_back(v);
                for (int c : dom_children[v]) {
                    if (!visited[c]) {
                        visited[c] = true;
                        bfs.push(c);
                    }
                }
            }
        }

        for (int i = (int)order.size() - 1; i >= 0; i--) {
            int v = order[i];

            if (v < n_real) {
                subtree_vv[v] = 1;
            } else {
                subtree_vv[v] = 0;
            }
            for (int c : dom_children[v]) {
                subtree_vv[v] += subtree_vv[c];
            }
        }

        for (int i = 0; i < n_virtual; i++) {
            int w = n_real + 1 + i;
            if (idom[w] == -1) continue;

            Edge e = virtual_to_edge[i];
            delta[e] += subtree_vv[w];
        }
    }
};

#endif
