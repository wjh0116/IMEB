#ifndef __EDGEBLOCK_H__
#define __EDGEBLOCK_H__

#include "infgraph.h"
#include "iheap.h"
#include "DESCE.h"
#include "CriticalEdge.h"

class EdgeBlock {
public:
    InfGraph& g;
    vector<int> seed_set;
    unsigned int k;

    vector<Edge> blocked_edge_set;
    double original_influence;
    double final_influence;
    double running_time;

    EdgeBlock(InfGraph& g, vector<int>& seeds, unsigned int k)
        : g(g), seed_set(seeds), k(k) {
        g.setSeedSet(seeds);
    }

    void greedyMC(int num_simulations = 10000) {
        auto start = chrono::high_resolution_clock::now();

        original_influence = g.estimateInfluence(num_simulations);
        cout << "Original influence: " << original_influence << endl;

        for (unsigned int iter = 0; iter < k; iter++) {
            int best_u = -1, best_i = -1;
            double best_inf = 1e18;

            for (unsigned int u = 0; u < g.n; u++) {
                for (size_t i = 0; i < g.adj[u].size(); i++) {
                    if (g.blocked[u][i]) continue;

                    g.blockByIdx(u, i);
                    double inf = g.estimateInfluence(num_simulations);
                    g.unblockByIdx(u, i);

                    if (inf < best_inf) {
                        best_inf = inf;
                        best_u = u;
                        best_i = i;
                    }
                }
            }

            if (best_u == -1) break;

            g.blockByIdx(best_u, best_i);
            blocked_edge_set.push_back({best_u, g.adj[best_u][best_i]});
            cout << "Blocked edge (" << best_u << ", " << g.adj[best_u][best_i]
                 << ") with influence after blocking: " << best_inf << endl;
        }

        final_influence = g.estimateInfluence(num_simulations);

        auto end = chrono::high_resolution_clock::now();
        running_time = chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000.0;
    }

    void greedyDT(int theta = 10000) {
        auto start = chrono::high_resolution_clock::now();

        original_influence = g.estimateInfluence(10000);
        cout << "Original influence (MC estimate): " << original_influence << endl;

        DESCE desce(g, seed_set, theta);
        for (unsigned int iter = 0; iter < k; iter++) {
            cout << "\n--- Greedy iteration " << (iter + 1) << "/" << k << " ---" << endl;

            desce.run();

            Edge best_edge = {-1, -1};
            double best_delta = -1.0;
            for (auto& [e, val] : desce.delta) {
                if (g.isBlocked(e.first, e.second)) continue;
                if (val > best_delta) {
                    best_delta = val;
                    best_edge = e;
                }
            }

            if (best_edge.first == -1) break;

            g.blockEdge(best_edge.first, best_edge.second);
            blocked_edge_set.push_back(best_edge);
            cout << "Blocked edge (" << best_edge.first << ", " << best_edge.second
                 << ") with estimated decrease: " << best_delta << endl;
        }

        final_influence = g.estimateInfluence(10000);

        auto end = chrono::high_resolution_clock::now();
        running_time = chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000.0;
    }

    void greedyReplaceDT(int theta = 10000) {
        auto start = chrono::high_resolution_clock::now();

        original_influence = g.estimateInfluence(10000);
        cout << "Original influence (MC estimate): " << original_influence << endl;

        set<int> in_seed(seed_set.begin(), seed_set.end());
        set<Edge> cand;
        for (int s : seed_set)
            for (size_t i = 0; i < g.adj[s].size(); i++)
                if (!in_seed.count(g.adj[s][i])) cand.insert({s, g.adj[s][i]});

        DESCE desce(g, seed_set, theta);
        size_t n1 = min((size_t)k, cand.size());
        for (size_t iter = 0; iter < n1; iter++) {
            cout << "\n--- Phase 1 iteration " << (iter + 1) << "/" << n1 << " ---" << endl;

            desce.run();

            Edge best_edge = {-1, -1};
            double best_delta = -1.0;
            for (const Edge& e : cand) {
                auto it = desce.delta.find(e);
                double val = (it == desce.delta.end()) ? 0.0 : it->second;
                if (val > best_delta) {
                    best_delta = val;
                    best_edge = e;
                }
            }

            if (best_edge.first == -1) break;

            g.blockEdge(best_edge.first, best_edge.second);
            blocked_edge_set.push_back(best_edge);
            cand.erase(best_edge);
            cout << "Blocked seed out-edge (" << best_edge.first << ", " << best_edge.second
                 << ") with estimated decrease: " << best_delta << endl;
        }

        vector<Edge> phase1 = blocked_edge_set;
        for (int idx = (int)phase1.size() - 1; idx >= 0; idx--) {
            Edge old_edge = phase1[idx];
            g.unblockEdge(old_edge.first, old_edge.second);

            desce.run();

            Edge best_edge = {-1, -1};
            double best_delta = -1.0;
            for (auto& [e, val] : desce.delta) {
                if (g.isBlocked(e.first, e.second)) continue;
                if (val > best_delta) {
                    best_delta = val;
                    best_edge = e;
                }
            }

            auto it = desce.delta.find(old_edge);
            double old_delta = (it == desce.delta.end()) ? 0.0 : it->second;
            if (best_edge.first == -1 || old_delta >= best_delta) best_edge = old_edge;

            g.blockEdge(best_edge.first, best_edge.second);
            if (best_edge == old_edge) {
                cout << "Phase 2: kept (" << old_edge.first << ", " << old_edge.second << ")" << endl;
                break;
            }

            *find(blocked_edge_set.begin(), blocked_edge_set.end(), old_edge) = best_edge;
            cout << "Phase 2: replaced (" << old_edge.first << ", " << old_edge.second
                 << ") with (" << best_edge.first << ", " << best_edge.second << ")" << endl;
        }

        final_influence = g.estimateInfluence(10000);

        auto end = chrono::high_resolution_clock::now();
        running_time = chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000.0;
    }

    void greedyDTplus(int theta = 10000) {
        auto start = chrono::high_resolution_clock::now();

        original_influence = g.estimateInfluence(10000);
        cout << "Original influence (MC estimate): " << original_influence << endl;

        CriticalEdge ce(g, seed_set, theta);
        for (unsigned int iter = 0; iter < k; iter++) {
            cout << "\n--- Greedy iteration " << (iter + 1) << "/" << k << " ---" << endl;

            ce.run();

            Edge best_edge = {-1, -1};
            long long best_count = 0;
            for (auto& [e, cnt] : ce.count) {
                if (cnt > best_count) {
                    best_count = cnt;
                    best_edge = e;
                }
            }

            if (best_edge.first == -1) break;

            g.blockEdge(best_edge.first, best_edge.second);
            blocked_edge_set.push_back(best_edge);
            cout << "Blocked edge (" << best_edge.first << ", " << best_edge.second
                 << ") with weighted coverage: " << best_count << endl;
        }

        final_influence = g.estimateInfluence(10000);

        auto end = chrono::high_resolution_clock::now();
        running_time = chrono::duration_cast<chrono::milliseconds>(end - start).count() / 1000.0;
    }

    void heuristicProbDeg() {
        auto start = chrono::high_resolution_clock::now();
        original_influence = g.estimateInfluence(10000);
        cout << "Original influence (MC estimate): " << original_influence << endl;

        vector<pair<double, pair<int,int>>> cand;
        cand.reserve(g.m);
        for (unsigned int u = 0; u < g.n; u++)
            for (size_t i = 0; i < g.adj[u].size(); i++) {
                int v = g.adj[u][i];
                cand.push_back({g.prob[u][i] * (double)g.out_deg[v], {(int)u, (int)i}});
            }
        blockTopK(cand);

        final_influence = g.estimateInfluence(10000);
        running_time = chrono::duration_cast<chrono::milliseconds>(
                           chrono::high_resolution_clock::now() - start).count() / 1000.0;
    }

    void heuristicSeedAdj() {
        auto start = chrono::high_resolution_clock::now();
        original_influence = g.estimateInfluence(10000);
        cout << "Original influence (MC estimate): " << original_influence << endl;

        vector<pair<double, pair<int,int>>> cand;
        for (int s : seed_set)
            for (size_t i = 0; i < g.adj[s].size(); i++)
                cand.push_back({g.prob[s][i], {s, (int)i}});
        blockTopK(cand);

        final_influence = g.estimateInfluence(10000);
        running_time = chrono::duration_cast<chrono::milliseconds>(
                           chrono::high_resolution_clock::now() - start).count() / 1000.0;
    }

    void blockTopK(vector<pair<double, pair<int,int>>>& cand) {
        size_t kk = min((size_t)k, cand.size());
        partial_sort(cand.begin(), cand.begin() + kk, cand.end(),
                     [](const pair<double,pair<int,int>>& a,
                        const pair<double,pair<int,int>>& b) { return a.first > b.first; });
        for (size_t j = 0; j < kk; j++) {
            int u = cand[j].second.first, i = cand[j].second.second;
            g.blockByIdx(u, i);
            blocked_edge_set.push_back({u, g.adj[u][i]});
        }
    }

    void printResults() {
        cout << "======================================" << endl;
        cout << "Results:" << endl;
        cout << "Algorithm: " << args.algo << endl;
        cout << "k (edges blocked): " << k << endl;
        cout << "Original influence: " << original_influence << endl;
        cout << "Final influence: " << final_influence << endl;
        cout << "Influence reduction: " << (original_influence - final_influence) << endl;
        cout << "Reduction ratio: " << (original_influence - final_influence) / original_influence * 100 << "%" << endl;
        cout << "Running time: " << running_time << " seconds" << endl;
        cout << "Blocked edges:" << endl;
        for (const auto& e : blocked_edge_set) {
            cout << "  (" << e.first << ", " << e.second << ")" << endl;
        }
        cout << "======================================" << endl;
    }

    void saveResults(string filename) {
        ofstream fout(filename);
        fout << "algorithm=" << args.algo << endl;
        fout << "k=" << k << endl;
        fout << "original_influence=" << original_influence << endl;
        fout << "final_influence=" << final_influence << endl;
        fout << "influence_reduction=" << (original_influence - final_influence) << endl;
        fout << "running_time=" << running_time << endl;
        fout << "blocked_edges:" << endl;
        for (const auto& e : blocked_edge_set) {
            fout << e.first << " " << e.second << endl;
        }
        fout.close();
    }
};

#endif
