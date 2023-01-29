#include <iostream>
#include <vector>
#include <set>
#include <queue>
using namespace std;

const int INF = 1'000'000'000;
const double CALC_TIME = 5.0;

vector<int> U, V, W;

struct edge {
    int to, cost;
    edge(int to, int cost) : to(to), cost(cost) {}
};

using graph = vector<set<int>>;
graph& operator +=(graph &self, const set<int> &S) {
    for (auto i : S) {
        int u = U[i], v = V[i];
        self[u].insert(v);
        self[v].insert(u);
    }
    return self;
};
graph& operator -=(graph &self, const set<int> &S) {
    for (auto i : S) {
        int u = U[i], v = V[i];
        self[u].erase(v);
        self[v].erase(u);
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

    auto base_ms = 0;
    auto tic = [&base_ms]() { base_ms = clock(); };
    auto toc = [&base_ms]() { return clock() - base_ms; };

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

    auto calc = [&](const vector<set<int>> &S) -> long long {
        long long res = 0;

        for (int d = 0; d < (int)S.size(); ++d) {
            G -= S[d];

            int s = rand() % N;
            auto ndist = dijkstra(s, G);
                
            for (int g = 0; g < N; ++g) {
                res += ndist[g] - dist[s][g];
            }

            G += S[d];
        }

        return res;
    };

    vector best(D, set<int>());
    for (int i = 0; i < M; ++i) {
        best[i % D].insert(i);
    }

    long long best_min = calc(best);

    int loop_count = 0;

    tic();
    while (toc() < CALC_TIME * CLOCKS_PER_SEC) {
        vector S(D, set<int>());
        for (int i = 0; i < M; ++i) {
            int d = rand() % D;
            while (S[d].size() >= K) d = rand() % D;
            S[d].insert(i);
        }

        long long res = calc(S);
        if (res < best_min) {
            best_min = res;
            best = S;
        }

        ++loop_count;
    }

    output(M, best);
    cerr << "loop = " << loop_count << endl;

    return 0;
}
