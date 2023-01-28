#include <iostream>
#include <vector>
#include <set>
#include <queue>
using namespace std;

const int INF = 1'000'000'000;

struct edge {
    int to, cost;
    edge(int to, int cost) : to(to), cost(cost) {}
};

vector<int> dijkstra(const int s, const vector<vector<edge>> &G) {
    int N = G.size();

    vector dist(N, INF);
    dist[s] = 0;

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> q;
    q.push({0, s});

    while (q.size()) {
        auto [d, v] = q.top(); q.pop();
        if (dist[v] < d) continue;
        for (auto [to, cost] : G[v]) {
            if (dist[v] + cost < dist[to]) {
                dist[to] = dist[v] + cost;
                q.push({dist[to], to});
            }
        }
    }

    return dist;
}

long long calc(const vector<set<int>> &S) {
    
}

int main() {
    int N, M, D, K;
    cin >> N >> M >> D >> K;

    vector<int> U(M), V(M), W(M);
    for (int i = 0; i < M; ++i) {
        cin >> U[i] >> V[i] >> W[i];
        --U[i], --V[i];
    }

    vector G(N, vector<edge>());
    for (int i = 0; i < M; ++i) {
        int u = U[i], v = V[i], w = W[i];
        G[u].push_back(edge(v, w));
        G[v].push_back(edge(u, w));
    }

/*
    vector<int> X(N), Y(N);
    for (int i = 0; i < N; ++i) {
        cin >> X[i] >> Y[i];
    }
*/

    vector S(D, set<int>());
    for (int i = 0; i < M; ++i) {
        S[i % D].insert(i);
    }

    vector dist(N, vector<int>());
    for (int i = 0; i < N; ++i) {
        dist[i] = dijkstra(i, G);
    }

    long long res = 0;
    for (int d = 0; d < D; ++d) {
        vector nG(N, vector<edge>());
        for (int i = 0; i < M; ++i) {
            if (S[d].count(i)) continue;
            int u = U[i], v = V[i], w = W[i];
            G[u].push_back(edge(v, w));
            G[v].push_back(edge(u, w));
        }

        for (int i = 0; i < N; ++i) {
            auto ndist = dijkstra(i, nG);
            
            for (int j = i + 1; j < N; ++j) {
                res += ndist[j] - dist[i][j];
            }
        }
    }

    vector<int> R(M);
    for (int d = 0; d < D; ++d) {
        for (auto &i : S[d]) {
            R[i] = d;
        }
    }
    for (auto &r : R) cout << r + 1 << " ";

    cerr << "score = " << res << endl;

    return 0;
}
