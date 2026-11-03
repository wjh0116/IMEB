#ifndef __GRAPH_H__
#define __GRAPH_H__

#include "head.h"

class Graph {
public:
    unsigned int n, m;

    vector<int> out_deg;
    vector<int> in_deg;

    vector<vector<int>> adj;
    vector<vector<int>> rev_adj;

    vector<vector<double>> prob;
    vector<vector<double>> rev_prob;

    string folder;
    string graph_file;

    vector<vector<bool>> blocked;

    Graph(string folder, string graph_file) : folder(folder), graph_file(graph_file) {
        readAttribute();
        readGraph();
    }

    void readAttribute() {
        string attr_file = folder + "/attribute.txt";
        ifstream fin(attr_file);
        if (!fin.is_open()) {
            cerr << "Error: cannot open " << attr_file << endl;
            exit(1);
        }
        string line;
        while (getline(fin, line)) {
            if (line.substr(0, 2) == "n=") {
                n = stoi(line.substr(2));
            } else if (line.substr(0, 2) == "m=") {
                m = stoi(line.substr(2));
            }
        }
        fin.close();
        cout << "Graph: n=" << n << ", m=" << m << endl;
    }

    void readGraph() {
        adj.resize(n);
        rev_adj.resize(n);
        prob.resize(n);
        rev_prob.resize(n);
        out_deg.resize(n, 0);
        in_deg.resize(n, 0);

        string gfile = folder + "/" + graph_file;
        ifstream fin(gfile);
        if (!fin.is_open()) {
            cerr << "Error: cannot open " << gfile << endl;
            exit(1);
        }

        string first_line;
        getline(fin, first_line);
        istringstream iss(first_line);
        int u, v;
        double p;
        unsigned int edge_count = 0;

        if (iss >> u >> v >> p) {

            adj[u].push_back(v);
            prob[u].push_back(p);
            rev_adj[v].push_back(u);
            rev_prob[v].push_back(p);
            out_deg[u]++;
            in_deg[v]++;
            edge_count++;
        }

        while (fin >> u >> v >> p) {
            adj[u].push_back(v);
            prob[u].push_back(p);
            rev_adj[v].push_back(u);
            rev_prob[v].push_back(p);
            out_deg[u]++;
            in_deg[v]++;
            edge_count++;
        }
        fin.close();

        cout << "Loaded " << edge_count << " edges from " << gfile << endl;
        assert(edge_count == m);

        blocked.resize(n);
        for (unsigned int i = 0; i < n; i++) {
            blocked[i].assign(adj[i].size(), false);
        }
    }

    vector<int> readSeedSet(int rumor_num) {
        vector<int> seeds;
        string seed_file = folder + "/rumorSet_" + to_string(rumor_num) + ".txt";
        ifstream fin(seed_file);
        if (!fin.is_open()) {
            cerr << "Error: cannot open " << seed_file << endl;
            exit(1);
        }
        int node;
        while (fin >> node) {
            seeds.push_back(node);
        }
        fin.close();
        cout << "Loaded " << seeds.size() << " seed nodes from " << seed_file << endl;
        return seeds;
    }

    void blockByIdx(int u, int i) {
        blocked[u][i] = true;
    }

    void unblockByIdx(int u, int i) {
        blocked[u][i] = false;
    }

    bool blockEdge(int u, int v) {
        for (size_t i = 0; i < adj[u].size(); i++) {
            if (adj[u][i] == v && !blocked[u][i]) {
                blocked[u][i] = true;
                return true;
            }
        }
        return false;
    }

    bool unblockEdge(int u, int v) {
        for (size_t i = 0; i < adj[u].size(); i++) {
            if (adj[u][i] == v && blocked[u][i]) {
                blocked[u][i] = false;
                return true;
            }
        }
        return false;
    }

    bool isBlocked(int u, int v) const {
        for (size_t i = 0; i < adj[u].size(); i++) {
            if (adj[u][i] == v) return blocked[u][i];
        }
        return false;
    }

    double getEffectiveProb(int u, int idx) const {
        if (blocked[u][idx]) return 0.0;
        return prob[u][idx];
    }

    void resetBlocked() {
        for (unsigned int i = 0; i < n; i++) {
            fill(blocked[i].begin(), blocked[i].end(), false);
        }
    }

};

#endif
