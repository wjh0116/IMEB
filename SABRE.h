#ifndef __SABRE_H__
#define __SABRE_H__

#include "infgraph.h"
#include "ForwardSample.h"
#include "CriticalEdgeSet.h"
#include "AncestralEdge.h"
#include "Greedy.h"

class SABRE {
public:
    InfGraph& g;
    vector<int> seeds;
    int k;
    double eps, gamma, delta;
    double nminus;

    vector<Edge> Bstar;
    double rho_star = 0;
    double rhoL_hat = 0, rhoU_hat = 0, rhoUu_est = 0;
    double alpha_L = 0, alpha_U = 0;
    double r_lower = 0, tau_hat = 0;
    double rhoU_surrogate = 0;
    int    rounds_used = 0;
    long   theta_used = 0;
    long   theta_opt_fs = 0;
    long   theta_rho_used = 0;
    long   theta_max_v = 0, theta0_v = 0;
    int    imax_v = 0;
    double kappa_v = 0;
    double running_time = 0;
    double t_lower = 0, t_upper = 0, t_rho = 0;

    SABRE(InfGraph& g_, vector<int>& seeds_, int k_, double eps_, double gamma_, double delta_)
        : g(g_), seeds(seeds_), k(k_), eps(eps_), gamma(gamma_), delta(delta_),
          nminus((double)g_.n - (double)seeds_.size()),
          is_seed(g_.n, 0), gen(), fs(g_, seeds_, gen), ceb(g_),
          QL1(g_.n), QL2(g_.n), QU1(g_.n), QU2(g_.n) {
        gen.seed((uint32_t)chrono::system_clock::now().time_since_epoch().count());
        for (int s : seeds) is_seed[s] = 1;
    }

    void run() {
        auto t0 = chrono::high_resolution_clock::now();
        const double thr = 1.0 - 1.0 / M_E - eps;

        kappa_v = computeKappa();
        double a = 1.0 - 1.0 / M_E;
        double lnCmk = lnBinom((double)g.m, k);
        double ln12d = log(12.0 / delta);
        double num = a * sqrt(ln12d) + sqrt(a * (lnCmk + ln12d));
        double theta_max = 2.0 * nminus * num * num / (eps * eps * max(kappa_v, 1e-12));
        long theta0 = (long)llround(theta_max * eps * eps * max(kappa_v, 1e-12) / nminus);
        theta0 = max(1L, theta0);
        long tmax = max((long)ceil(theta_max), theta0);
        int imax = max(1, (int)ceil(log2((double)tmax / (double)theta0))) + 1;
        double zeta = log(12.0 * imax / delta);
        theta_max_v = tmax; theta0_v = theta0; imax_v = imax;

        long cur = theta0;
        drawInto(QL1, QU1, cur);
        drawInto(QL2, QU2, cur);

        GreedyResult gL, gU;
        for (int i = 1; i <= imax; i++) {
            rounds_used = i;
            if (alpha_L < thr) { auto tg = chrono::high_resolution_clock::now();
                gL = weightedGreedy(QL1, k); alpha_L = certAlpha(QL1, QL2, gL, zeta); tms_gl += now_ms(tg); }
            if (alpha_U < thr) { auto tg = chrono::high_resolution_clock::now();
                gU = weightedGreedy(QU1, k); alpha_U = certAlpha(QU1, QU2, gU, zeta); tms_gu += now_ms(tg); }
            if (alpha_L >= thr && alpha_U >= thr) break;

            bool growL = (alpha_L < thr) && (long)QL1.num_samples < tmax;
            bool growU = (alpha_U < thr) && (long)QU1.num_samples < tmax;
            if (!growL && !growU) break;
            if (growL && growU) {
                drawInto(QL1, QU1, cur); drawInto(QL2, QU2, cur);
            } else if (growL) {
                drawSide(QL1, true, cur);  drawSide(QL2, true, cur);
            } else {
                drawSide(QU1, false, cur); drawSide(QU2, false, cur);
            }
            cur *= 2;
        }
        theta_used = max((long)QL1.num_samples, (long)QU1.num_samples);

        auto trho = chrono::high_resolution_clock::now();
        estimateRhoBoth(gL.B, gU.B, rhoL_hat, rhoU_hat);
        if (rhoU_hat >= rhoL_hat) { Bstar = gU.B; rho_star = rhoU_hat; }
        else                      { Bstar = gL.B; rho_star = rhoL_hat; }

        rhoU_surrogate = rhoUu_est;
        double mfac = (1.0 - gamma) / (1.0 + gamma);
        tau_hat = (rhoU_surrogate > 0) ? rhoU_hat / rhoU_surrogate : 0.0;
        double tau_c = min(1.0, tau_hat);
        r_lower = (1.0 - 1.0 / M_E - eps) * mfac * mfac * tau_c;

        tms_rho = now_ms(trho);

        t_lower = (0.5 * tms_sample + tms_ce + tms_gl) / 1000.0;
        t_upper = (0.5 * tms_sample + tms_ae + tms_gu) / 1000.0;
        t_rho   = tms_rho / 1000.0;

        running_time = chrono::duration_cast<chrono::milliseconds>(
                           chrono::high_resolution_clock::now() - t0).count() / 1000.0;
    }

private:
    vector<char> is_seed;
    SFMTRng gen;
    ForwardSampler fs;
    CESetBuilder ceb;
    Pool QL1, QL2, QU1, QU2;
    vector<Edge> R, P;
    double tms_sample = 0, tms_ce = 0, tms_ae = 0, tms_gl = 0, tms_gu = 0, tms_rho = 0;
    static inline double now_ms(chrono::high_resolution_clock::time_point t) {
        return chrono::duration_cast<chrono::nanoseconds>(
                   chrono::high_resolution_clock::now() - t).count() / 1e6;
    }

    static double lnBinom(double m, int k) {
        return lgamma(m + 1.0) - lgamma(k + 1.0) - lgamma(m - k + 1.0);
    }

    double computeKappa() {
        vector<char> cand(g.n, 0);
        vector<int> cands;
        for (int s : seeds)
            for (int w : g.adj[s])
                if (!is_seed[w] && !cand[w]) { cand[w] = 1; cands.push_back(w); }

        vector<double> qv;
        qv.reserve(cands.size());
        for (int v : cands) {
            const auto& in = g.rev_adj[v];
            const auto& ip = g.rev_prob[v];
            int d = (int)in.size();
            double prodNZ = 1.0; int zeros = 0;
            for (int j = 0; j < d; j++) {
                double om = 1.0 - ip[j];
                if (om <= 0.0) zeros++; else prodNZ *= om;
            }
            double best = 0.0;
            for (int j = 0; j < d; j++) {
                if (!is_seed[in[j]]) continue;
                double om = 1.0 - ip[j];
                double prodExcl;
                if (om <= 0.0) prodExcl = (zeros == 1) ? prodNZ : 0.0;
                else           prodExcl = (zeros >= 1) ? 0.0 : prodNZ / om;
                best = max(best, ip[j] * prodExcl);
            }
            if (best > 0) qv.push_back(best);
        }
        sort(qv.rbegin(), qv.rend());
        double kap = 0.0;
        for (int i = 0; i < (int)qv.size() && i < k; i++) kap += qv[i];
        return kap;
    }

    void drawInto(Pool& CEpool, Pool& AEpool, long count) {
        long got = 0;
        while (got < count) {
            auto ts = chrono::high_resolution_clock::now();
            bool ok = fs.sample(); theta_opt_fs++;
            tms_sample += now_ms(ts);
            if (!ok) {
                R.clear(); CEpool.addSample(R, 0.0);
                P.clear(); AEpool.addSample(P, 0.0);
                got++;
                continue;
            }
            double N = fs.N;
            auto tc = chrono::high_resolution_clock::now();
            R.clear(); ceb.build(fs, R); CEpool.addSample(R, N);
            tms_ce += now_ms(tc);
            auto ta = chrono::high_resolution_clock::now();
            P.clear(); buildAESet(fs, P); AEpool.addSample(P, N);
            tms_ae += now_ms(ta);
            got++;
        }
    }

    void drawSide(Pool& pool, bool isCE, long count) {
        long got = 0;
        while (got < count) {
            auto t0 = chrono::high_resolution_clock::now();
            bool ok = fs.sample(); theta_opt_fs++;
            if (!ok) {
                if (isCE) { R.clear(); pool.addSample(R, 0.0); }
                else      { P.clear(); pool.addSample(P, 0.0); }
            } else if (isCE) {
                R.clear(); ceb.build(fs, R); pool.addSample(R, fs.N);
            } else {
                P.clear(); buildAESet(fs, P); pool.addSample(P, fs.N);
            }
            double dt = now_ms(t0);
            if (isCE) tms_ce += dt; else tms_ae += dt;
            got++;
        }
    }

    double certAlpha(Pool& Q1, Pool& Q2, const GreedyResult& gr, double zeta) {
        double t1 = Q1.num_samples, t2 = Q2.num_samples;
        double Lam2 = evalCoverage(Q2, gr.B) / nminus;
        double Lamu = gr.Lambda_u / nminus;
        double base = sqrt(Lam2 + 2.0 * zeta / 9.0) - sqrt(zeta / 2.0);
        double rho_l = (nminus / t2) * (base * base - zeta / 18.0);
        if (rho_l < 0) rho_l = 0;
        double su = sqrt(Lamu + zeta / 2.0) + sqrt(zeta / 2.0);
        double rho_u = (nminus / t1) * su * su;
        return rho_u > 0 ? rho_l / rho_u : 0.0;
    }

    void estimateRhoBoth(const vector<Edge>& BL, const vector<Edge>& BU,
                         double& rhoL, double& rhoU) {
        long long nn = g.n;
        auto enc = [nn](const Edge& e) { return (long long)e.first * nn + e.second; };
        unordered_set<long long> SL, SU;
        for (auto& e : BL) SL.insert(enc(e));
        for (auto& e : BU) SU.insert(enc(e));

        double dp = delta / 6.0;

        double Gamma = 2.0 * nminus * (1.0 + gamma) * (1.0 + gamma / 3.0) * log(2.0 / dp) / (gamma * gamma);

        double sumL = 0, sumU = 0, sumUu = 0;
        double thetaL = -1, thetaU = -1, thetaUu = -1;
        long theta = 0;
        vector<char> vis(g.n, 0); vector<int> bq, tch;
        vector<Edge> Pu;
        while (thetaL < 0 || thetaU < 0 || thetaUu < 0) {
            theta++;
            if (theta > (long)200 * theta_max_v + 1000000) break;
            if (!fs.sample()) continue;
            double zL = (thetaL < 0) ? zScore(SL, enc, vis, bq, tch) : 0.0;
            double zU = (thetaU < 0) ? zScore(SU, enc, vis, bq, tch) : 0.0;
            if (thetaL < 0) { sumL += zL; if (sumL >= Gamma) thetaL = theta - (sumL - Gamma) / max(zL, 1e-12); }
            if (thetaU < 0) { sumU += zU; if (sumU >= Gamma) thetaU = theta - (sumU - Gamma) / max(zU, 1e-12); }
            if (thetaUu < 0) {
                Pu.clear(); buildAESet(fs, Pu);
                bool hit = false;
                for (auto& e : Pu) if (SU.count(enc(e))) { hit = true; break; }
                double zUu = hit ? (double)fs.N : 0.0;
                sumUu += zUu; if (sumUu >= Gamma) thetaUu = theta - (sumUu - Gamma) / max(zUu, 1e-12);
            }
        }
        theta_rho_used = theta;
        rhoL = (thetaL > 0) ? Gamma / thetaL : 0.0;
        rhoU = (thetaU > 0) ? Gamma / thetaU : 0.0;
        rhoUu_est = (thetaUu > 0) ? Gamma / thetaUu : 0.0;
    }

    template <class Enc>
    double zScore(unordered_set<long long>& B, Enc enc,
                  vector<char>& vis, vector<int>& bq, vector<int>& tch) {

        for (int x : tch) vis[x] = 0;
        tch.clear();
        bq.clear();
        for (int s : seeds) if (!vis[s]) { vis[s] = 1; tch.push_back(s); bq.push_back(s); }
        bool reached = vis[fs.v];
        for (size_t h = 0; h < bq.size() && !reached; h++) {
            int u = bq[h];
            for (int w : fs.sadj[u]) {
                if (vis[w]) continue;
                if (B.count((long long)u * g.n + w)) continue;
                vis[w] = 1; tch.push_back(w); bq.push_back(w);
                if (w == fs.v) { reached = true; break; }
            }
        }
        return reached ? 0.0 : (double)fs.N;
    }
};

#endif
