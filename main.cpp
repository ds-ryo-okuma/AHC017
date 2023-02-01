#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <numeric>
#include <cmath>
#include <climits>
using namespace std;

const int INF = 1'000'000'000;
const double INITIAL_TIME = 1.0;
const double CALC_TIME = 4.5;

vector<int> U, V, W;

struct edge {
    int to, cost;
    edge(int to, int cost) : to(to), cost(cost) {}
};

using graph = vector<set<int>>;
graph& operator +=(graph &self, const set<int> &S) {
    for (auto i : S) {
        int u = U[i], v = V[i];
        self[u].insert(i);
        self[v].insert(i);
    }
    return self;
};
graph& operator -=(graph &self, const set<int> &S) {
    for (auto i : S) {
        int u = U[i], v = V[i];
        self[u].erase(i);
        self[v].erase(i);
    }
    return self;
};

vector<int> dijkstra(const int s, const graph &G) {
    int N = G.size();

    vector dist(N, INF);
    dist[s] = 0;

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> q;
    q.push({0, s});

    while (q.size()) {
        auto [d, v] = q.top(); q.pop();
        if (dist[v] < d) continue;
        for (auto i : G[v]) {
            int to = U[i] == v ? V[i] : U[i];
            int cost = W[i];
            if (dist[v] + cost < dist[to]) {
                dist[to] = dist[v] + cost;
                q.push({dist[to], to});
            }
        }
    }

    return dist;
}


unsigned int randxor() {
    static unsigned int x = 123456789, y = 362436069, z = 521288629, w = 88675123;
    unsigned int t;
    t = (x ^ (x << 11)); x = y; y = z; z = w; return(w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
}

static double rand01()
{
    return (randxor() + 0.5) * (1.0 / UINT_MAX);
}

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


long long comp(const vector<long long> &A, const vector<long long> &B) {
    return accumulate(A.begin(), A.end(), 0ll) - accumulate(B.begin(), B.end(), 0LL);
}


void output (int M, const vector<set<int>> &S) {
    vector<int> R(M);
    for (int d = 0; d < (int)S.size(); ++d) {
        for (auto &i : S[d]) {
            R[i] = d;
        }
    }
    for (auto &r : R) cout << r + 1 << " ";
}


int main() {
    int N, M, D, K;
    cin >> N >> M >> D >> K;

    U.resize(M), V.resize(M), W.resize(M);
    for (int i = 0; i < M; ++i) {
        cin >> U[i] >> V[i] >> W[i];
        --U[i], --V[i];
    }

    clock_t base_ms = 0;
    auto tic = [&base_ms]() { base_ms = clock(); };
    auto toc = [&base_ms]() -> clock_t { return clock() - base_ms; };

    graph G(N);
    for (int i = 0; i < M; ++i) {
        int u = U[i], v = V[i], w = W[i];
        G[u].insert(i);
        G[v].insert(i);
    }

    vector dist(N, vector<int>());
    for (int i = 0; i < N; ++i) {
        dist[i] = dijkstra(i, G);
    }

/*
    vector<int> X(N), Y(N);
    for (int i = 0; i < N; ++i) {
        cin >> X[i] >> Y[i];
    }
*/

    auto calc = [&](
        const vector<set<int>> &S,
        const vector<set<int>> &old_state = vector<set<int>>(),
        const vector<long long> &old_score = vector<long long>()
        ) -> vector<long long> {

        bool initial = old_state.empty();
        vector<long long> res;
        if (initial) res.resize(D);
        else res = old_score;

        for (int d = 0; d < D; ++d) {
            if (!initial && old_state[d] == S[d]) continue;

            res[d] = 0;

            G -= S[d];

            int s = randxor() % N;
            auto ndist = dijkstra(s, G);

            G += S[d];
                
            for (int g = 0; g < N; ++g) res[d] += ndist[g] - dist[s][g];
        }

        return res;
    };


    vector best(D, set<int>());
    vector<long long> score(D, 1'000'000'000'000'000ll);

    vector S(D, set<int>());
    for (int i = 0; i < M; ++i) {
        S[i % D].insert(i);
    }

    vector<long long> res = calc(S);
    best = S; score = res;

    int initial_count = 0;

    while (toc() < INITIAL_TIME * CLOCKS_PER_SEC) {
        S = make_random(D, M, K);

        res = calc(S);
        if (comp(res, score) < 0) {
            score = res;
            best = S;
        }

        ++initial_count;
    }

    int loop_count = 0;

    const static double START_TEMP = 1'000'000; // 開始時の温度
    const static double END_TEMP   =  100;  // 終了時の温度
    const static double END_TIME   =  CALC_TIME * CLOCKS_PER_SEC; // 終了時間（秒）

    tic();
    
    clock_t cur = 0;
    while ((cur = toc()) < CALC_TIME * CLOCKS_PER_SEC) {
        const double progressRatio = cur / END_TIME;
        const double temp = START_TEMP + (END_TEMP - START_TEMP) * progressRatio;

        S = neighbor(K, best);

        res = calc(S, best, score);

        long long diff = comp(score, res);
        const double probability = exp(diff / temp);
        bool FORCE_NEXT = probability > rand01();

        if (FORCE_NEXT) {
            score = res;
            best = S;
        }

        ++loop_count;
    }

    output(M, best);
    cerr << "INITIAL_COUNT = " << initial_count << endl;
    cerr << "LOOP_COUNT = " << loop_count << endl;

    return 0;
}
