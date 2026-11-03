#ifndef __HEAD_H__
#define __HEAD_H__

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <numeric>
#include <random>
#include <chrono>

#include "sfmt/SFMT.h"

using namespace std;

class SFMTRng {
public:
    sfmt_t sfmt;
    SFMTRng() { sfmt_init_gen_rand(&sfmt, 0); }
    void seed(uint32_t s) { sfmt_init_gen_rand(&sfmt, s); }

    double randomDouble() { return sfmt_genrand_real2(&sfmt); }

    uint32_t randomUint32() { return sfmt_genrand_uint32(&sfmt); }
};

string __head_version;

typedef pair<int, int> Edge;

class Argument {
public:
    unsigned int k;
    string dataset;
    string res;
    string algo;
    int rumor_num;
    double delta;
    double gamma;
    double epsilon;
    int theta;
};

Argument args;

void parseArgs(int argn, char** argv) {
    for (int i = 1; i < argn; i++) {
        if (string(argv[i]) == "-dataset") args.dataset = string(argv[++i]);
        else if (string(argv[i]) == "-k") args.k = atoi(argv[++i]);
        else if (string(argv[i]) == "-rumorNum") args.rumor_num = atoi(argv[++i]);
        else if (string(argv[i]) == "-algo") args.algo = string(argv[++i]);
        else if (string(argv[i]) == "-epsilon") args.epsilon = atof(argv[++i]);
        else if (string(argv[i]) == "-gamma") args.gamma = atof(argv[++i]);
        else if (string(argv[i]) == "-delta") args.delta = atof(argv[++i]);
        else if (string(argv[i]) == "-theta") args.theta = atoi(argv[++i]);
        else if (string(argv[i]) == "-res") args.res = string(argv[++i]);
    }
}

class OutputInfo {
public:
    OutputInfo(int argn, char** argv) {
        cout << "======================================" << endl;
        cout << "Influence Minimization via Edge Blocking" << endl;
        cout << "Version: " << __head_version << endl;
        cout << "Arguments: ";
        for (int i = 0; i < argn; i++) cout << argv[i] << " ";
        cout << endl;
        cout << "======================================" << endl;
    }
};

#endif
