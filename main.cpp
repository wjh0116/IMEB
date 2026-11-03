#include "EdgeBlock.h"
#include "SABRE.h"

void Run(int argn, char** argv) {
    parseArgs(argn, argv);

    cout << "Dataset: " << args.dataset << endl;
    cout << "Algorithm: " << args.algo << endl;
    cout << "k: " << args.k << endl;
    cout << "Epsilon: " << args.epsilon << endl;

    InfGraph g(args.dataset, "graph.txt");
    vector<int> seeds = g.readSeedSet(args.rumor_num);
    g.setSeedSet(seeds);

    if (args.algo == "SABRE") {
        double eps   = args.epsilon > 0 ? args.epsilon : 0.1;
        double gamma = args.gamma   > 0 ? args.gamma   : 0.3;
        double delta = args.delta   > 0 ? args.delta   : 1.0 / g.n;
        SABRE sabre(g, seeds, args.k, eps, gamma, delta);
        sabre.run();

        double sigmaS = g.estimateInfluence(10000);
        for (auto& e : sabre.Bstar) g.blockEdge(e.first, e.second);
        double sigmaB = g.estimateInfluence(10000);
        g.resetBlocked();

        cout << "======================================" << endl;
        cout << "Algorithm: " << args.algo << endl;
        cout << "k (edges blocked): " << args.k << endl;
        cout << "Original influence: " << sigmaS << endl;
        cout << "Final influence: " << sigmaB << endl;
        cout << "Influence reduction: " << (sigmaS - sigmaB) << endl;
        cout << "Reduction ratio: " << (sigmaS - sigmaB) / sigmaS * 100 << "%" << endl;
        cout << "Approximation guarantee (Thm): r_lower=" << sabre.r_lower
             << " tau=" << sabre.tau_hat << endl;
        cout << "(internal: alpha_L=" << sabre.alpha_L
             << " alpha_U=" << sabre.alpha_U << ")" << endl;
        cout << "rho_hat(B*): " << sabre.rho_star << endl;
        cout << "Runtime breakdown: lower-maximize=" << sabre.t_lower
             << "s upper-maximize=" << sabre.t_upper
             << "s rho-estimation=" << sabre.t_rho << "s" << endl;
        cout << "rounds=" << sabre.rounds_used << " theta_used=" << sabre.theta_used
             << " theta_max=" << sabre.theta_max_v << " kappa=" << sabre.kappa_v << endl;

        long opt_samples = sabre.theta_opt_fs;
        cout << "Samples (forward BFS): optimize=" << opt_samples
             << " (max per-pool theta=" << sabre.theta_used << ")"
             << " rho-est=" << sabre.theta_rho_used
             << " total=" << (opt_samples + sabre.theta_rho_used) << endl;
        cout << "Running time: " << sabre.running_time << " seconds" << endl;
        cout << "======================================" << endl;

        if (!args.res.empty()) {
            ofstream fout(args.res);

            fout << "algorithm=" << args.algo << endl;
            fout << "k=" << args.k << endl;
            fout << "original_influence=" << sigmaS << endl;
            fout << "final_influence=" << sigmaB << endl;
            fout << "influence_reduction=" << (sigmaS - sigmaB) << endl;
            fout << "running_time=" << sabre.running_time << endl;

            fout << "r_lower=" << sabre.r_lower << endl;
            fout << "tau_hat=" << sabre.tau_hat << endl;
            fout << "rho_star=" << sabre.rho_star << endl;
            fout << "alpha_L=" << sabre.alpha_L << endl;
            fout << "alpha_U=" << sabre.alpha_U << endl;
            fout << "t_lower=" << sabre.t_lower << endl;
            fout << "t_upper=" << sabre.t_upper << endl;
            fout << "t_rho=" << sabre.t_rho << endl;
            fout << "rounds=" << sabre.rounds_used << endl;
            fout << "theta_used=" << sabre.theta_used << endl;
            fout << "theta_max=" << sabre.theta_max_v << endl;
            fout << "kappa=" << sabre.kappa_v << endl;
            fout << "samples_optimize=" << sabre.theta_opt_fs << endl;
            fout << "samples_rho=" << sabre.theta_rho_used << endl;
            fout << "samples_total=" << (sabre.theta_opt_fs + sabre.theta_rho_used) << endl;
            fout << "blocked_edges:" << endl;
            for (auto& e : sabre.Bstar) fout << e.first << " " << e.second << endl;
            fout.close();
        }
        return;
    }

    EdgeBlock eb(g, seeds, args.k);

    if (args.algo == "Gre-MC") {
        eb.greedyMC(args.theta > 0 ? args.theta : 10000);
    } else if (args.algo == "Gre-DT") {
        eb.greedyDT(args.theta > 0 ? args.theta : 10000);
    } else if (args.algo == "Gre*-DT") {
        eb.greedyReplaceDT(args.theta > 0 ? args.theta : 10000);
    } else if (args.algo == "Gre-DT+") {
        eb.greedyDTplus(args.theta > 0 ? args.theta : 10000);
    } else if (args.algo == "ProbDeg") {
        eb.heuristicProbDeg();
    } else if (args.algo == "SeedAdj") {
        eb.heuristicSeedAdj();
    } else {
        cerr << "Unknown algorithm: " << args.algo << endl;
        cerr << "Available: SABRE, Gre-MC, Gre-DT, Gre*-DT, Gre-DT+, ProbDeg, SeedAdj" << endl;
        exit(1);
    }

    eb.printResults();

    if (!args.res.empty()) {
        eb.saveResults(args.res);
    }
}

int main(int argn, char** argv) {
    __head_version = "v0.1";
    OutputInfo info(argn, argv);
    Run(argn, argv);
    return 0;
}
