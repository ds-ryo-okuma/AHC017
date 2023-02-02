#include <vector>
#include <set>
#include<queue>

const int INF = 1'000'000'000;

class edge {
public:
    edge(int from, int to, int cost);
    int from, to, cost;
};
edge::edge(int from, int to, int cost) : from(from), to(to), cost(cost) {}



class graph {
public:
    int V, E;
    std::vector<std::set<int>> G;
    std::vector<edge> edges;

    graph(int V, int E) : V(V), E(E), G(V, std::set<int>()), edges() {};

    graph& operator +=(const std::set<int> &S);
    graph& operator -=(const std::set<int> &S);

    void add_edge(int from, int to, int cost);
    std::vector<int> dijkstra(const int s) const;
};

graph& graph::operator +=(const std::set<int> &S) {
    for (auto i : S) {
        auto [from, to, cost] = edges[i];
        G[from].insert(i);
        G[to].insert(i);
    }
    return *this;
}

graph& graph::operator -=(const std::set<int> &S) {
    for (auto i : S) {
        auto [from, to, cost] = edges[i];
        G[from].erase(i);
        G[to].erase(i);
    }
    return *this;
}

void graph::add_edge (int from, int to, int cost) {
    edges.emplace_back(from, to, cost);
    G[from].emplace(to);
    G[to].emplace(from);
}

std::vector<int> graph::dijkstra(const int s) const {
    std::vector dist(V, INF);
    dist[s] = 0;

    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> q;
    q.emplace(0, s);

    while (q.size()) {
        auto [d, v] = q.top(); q.pop();
        if (dist[v] < d) continue;
        for (auto i : G[v]) {
            auto [from, to, cost] = edges[i];
            if (v == to) std::swap(from, to);
            if (dist[v] + cost < dist[to]) {
                dist[to] = dist[v] + cost;
                q.emplace(dist[to], to);
            }
        }
    }

    return dist;
}
