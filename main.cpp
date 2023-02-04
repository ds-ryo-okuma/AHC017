#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <numeric>
#include <cmath>

#include "util.hpp"
#include "graph.hpp"
using namespace std;


const double INITIAL_TIME = 2.0;
const double CALC_TIME = 3.5;
const int CALC_NUM  = 1;

const static double START_TEMP = 1'000'000; // 開始時の温度
const static double END_TEMP   =  100;  // 終了時の温度
const static double END_TIME   =  CALC_TIME * CLOCKS_PER_SEC; // 終了時間（秒）

using state = set<int>;
using states = vector<state>;

using score_type = vector<double>;

vector<set<int>> make_random(int D, int M, int K) {
    vector res(D, set<int>());
    for (int i = 0; i < M; ++i) {
        int d = randxor() % D;
        while (res[d].size() >= K) d = randxor() % D;
        res[d].insert(i);
    }
    return res;
}


vector<set<int>> neighbor(int K, const vector<set<int>> &S) {
    int D = S.size();

    vector res(S);

    int i = randxor() % D, j = randxor() % (D - 1);
    if (i <= j) j++;

    // S_i と S_j を作り直す
    set<int> nS[2];
    for (auto s : res[i]) {
        int x = randxor() % 2;
        if (K <= nS[x].size()) x = 1 - x;
        nS[x].insert(s);
    }
    for (auto s : res[j]) {
        int x = randxor() % 2;
        if (K <= nS[x].size()) x = 1 - x;
        nS[x].insert(s);
    }

    res[i] = nS[0], res[j] = nS[1];

    return res;
}


score_type calc_score (const vector<set<int>> &S, graph &G, vector<vector<int>> &dist) {
    int D = S.size();
    score_type res(S.size());

    for (int d = 0; d < D; ++d) {
        G -= S[d];

        vector used(G.V, false);
        for (int x = 0; x < CALC_NUM; ++x) {
            int s = randxor() % G.V;
            while (used[s]) s = randxor() % G.V;
            used[s] = true;
            auto ndist = G.dijkstra(s);

            for (int g = 0; g < G.V; ++g) res[d] += (ndist[g] - dist[s][g]) / G.V;
        }

        G += S[d];
    }

    return res;
}


score_type calc_score (const vector<set<int>> &next, graph &G, vector<vector<int>> &dist, const vector<set<int>> &now, const score_type &score) {
    int D = next.size();
    score_type res(score);

    for (int d = 0; d < D; ++d) {
        if (next[d] == now[d]) continue;

        G -= next[d];

        vector used(G.V, false);
        for (int x = 0; x < CALC_NUM; ++x) {
            int s = randxor() % G.V;
            while (used[s]) s = randxor() % G.V;
            used[s] = true;
            auto ndist = G.dijkstra(s);

            for (int g = 0; g < G.V; ++g) res[d] += (ndist[g] - dist[s][g]) / G.V;
        }

        G += next[d];
    }

    return res;
}


int main() {
    int N, M, D, K;
    cin >> N >> M >> D >> K;

    graph all_path(N, M);

    for (int i = 0; i < M; ++i) {
        int u, v, w;
        cin >> u >> v >> w;
        all_path.add_edge(--u, --v, w);
    }

    vector dist(N, vector<int>());
    for (int i = 0; i < N; ++i) {
        dist[i] = all_path.dijkstra(i);
    }


    vector now(D, set<int>());
    for (int i = 0; i < M; ++i) {
        now[i % D].insert(i);
    }

    score_type score = calc_score(now, all_path, dist);

    int initial_count = 0;
    tic();
    while (toc() < INITIAL_TIME * CLOCKS_PER_SEC) {
        auto next = make_random(D, M, K);

        auto res = calc_score(next, all_path, dist);
        if (comp(res, score) < 0) {
            score = res;
            now = next;
        }

        ++initial_count;
    }


    int loop_count = 0;

    clock_t cur = 0;
    tic();
    while ((cur = toc()) < CALC_TIME * CLOCKS_PER_SEC) {
        const double progressRatio = cur / END_TIME;
        const double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

        auto next = neighbor(K, now);

        auto res = calc_score(next, all_path, dist, now, score);

        long long diff = comp(score, res);
        const double probability = exp(diff / temp);
        bool FORCE_NEXT = probability > rand01();

        if (FORCE_NEXT) {
            score = res;
            now = next;
        }

        ++loop_count;
    }

    output(M, now);
    cerr << "INITIAL_COUNT = " << initial_count << endl;
    cerr << "LOOP_COUNT = " << loop_count << endl;

    return 0;
}
