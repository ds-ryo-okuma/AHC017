#pragma once

#include <vector>
#include <set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <assert.h>

const int INF = 1'000'000'000;

class edge {
public:
    edge();
    edge(int from, int to, int cost);
    int from, to, cost;
};
edge::edge() {}
edge::edge(int from, int to, int cost) : from(from), to(to), cost(cost) {}



class graph {
    using pque = std::priority_queue<std::pair<long long, int>>;
public:
    // 頂点数
    int V;
    
    // 辺数
    int E;

    // グラフの隣接リスト G[i] := (隣接する頂点, 辺のコスト)
    std::vector<std::vector<std::pair<int, int>>> G;

    // 最短距離 dist[s][v] := 頂点sからvへの最短距離
    std::unordered_map<int, std::vector<long long>> dist, dist_bk;

    // 始点をsとするときの各頂点の最短距離木での親
    std::unordered_map<int, std::vector<int>> parent, parent_bk;

    // 追加・削除する辺のキュー
    std::queue<edge> add_edge_que, remove_edge_que;

    graph();
    graph(int V, int E);

    void add_edge(int from, int to, int cost);
    void dijkstra(int s);
    void dijkstra_sub_func(int s, pque &q);
    void dijkstra_add_edge(edge e);
    void dijkstra_remove_edge(edge e);
    void dfs(int s, int v, edge e, pque &q);
    void dfs_sub_tree(int s, int v, long long diff);
    long long sum();
    void backup();
    void commit();
    void rollback();
};

/**
 * @brief Construct a new graph::graph object
 * 
 */
graph::graph() {}

/**
 * @brief Construct a new graph::graph object
 * 
 * @param V 頂点数
 * @param E 辺数
 */
graph::graph(int V, int E) : V(V), G(V) {}

void graph::add_edge(int from, int to, int cost) {
    G[from].emplace_back(to, cost);
    G[to].emplace_back(from, cost);
}

/**
 * @brief 優先度付きキューが空になるまで dijkstra を行う
 * 
 * @param s 始点
 * @param q 優先度付きキュー
 */
void graph::dijkstra_sub_func(int s, pque &q) {
    auto &dist = this->dist[s];
    auto &parent = this->parent[s];

    while (q.size()) {
        auto [d, v] = q.top(); q.pop();
        d = -d;
        if (dist[v] < d) continue;

        for (auto [to, cost] : G[v]) {
            if (dist[v] + cost < dist[to]) {
                dist[to] = dist[v] + cost;
                parent[to] = v;
                q.emplace(-dist[to], to);
            }
        }
    }
}

/**
 * @brief 始点 s から各頂点への最短距離を求める
 * 
 * @param s 
 */
void graph::dijkstra(int s) {
    dist[s].assign(V, INF);
    dist[s][s] = 0;

    parent[s].assign(V, -1);

    pque q;
    q.emplace(0, s);

    dijkstra_sub_func(s, q);
}

/**
 * @brief 1 辺追加する場合の更新処理
 * 
 * @param e 追加する辺
 */
void graph::dijkstra_add_edge(edge e) {
    auto [v0, v1, cost] = e;

    add_edge(v0, v1, cost);
    add_edge_que.emplace(e);

    for (auto& [s, dist] : dist) {
        auto &parent = this->parent[s];

        pque q;

        if (dist[v0] + cost < dist[v1]) {
            // v0 -> v1 へ最短経路を更新
            long long old_dist = dist[v1], new_dist = dist[v0] + cost;
            
            // v1 の親を更新
            parent[v1] = v0;

            // 部分木を再帰的に更新
            dfs_sub_tree(s, v1, new_dist - old_dist);
        }
        if (dist[v1] + cost < dist[v0]) {
            // v1 -> v0 へ最短経路を更新
            long long old_dist = dist[v0], new_dist = dist[v1] + cost;
            
            // v0 の親を更新
            parent[v0] = v1;

            // 部分木を再帰的に更新
            dfs_sub_tree(s, v0, new_dist - old_dist);
        }
    }
}

/**
 * @brief 頂点vの部分木について、再計算の対象にする
 * 
 * @param s 始点
 * @param v 部分木の根
 * @param e 削除する辺
 * @param q dijkstraで使用する優先度付きキュー
 */
void graph::dfs(int s, int v, edge e, pque &q) {
    assert(v < V);

    auto &dist = this->dist[s];
    auto &parent = this->parent[s];

    // 頂点 v を未確定にする
    dist[v] = INF;
    parent[v] = -1;

    // 部分木に沿って再帰的に未確定にしていく
    for (auto [to, cost] : G[v]) {
        if (parent[to] != v) continue; 
        dfs(s, to, e, q);
    }

    // 頂点 v に隣接している頂点で、最短距離が確定している頂点をキューに追加
    for (auto [to, cost] : G[v]) {
        if (dist[to] < INF) {
            q.emplace(-dist[to], to);
        }
    }
}

/**
 * @brief 頂点 v の部分木を +diff で更新
 * 
 * @param v 
 * @param dist 
 * @param diff 
 * @param rooted_tree 
 */
void graph::dfs_sub_tree(int s, int v, long long diff) {
    auto &dist = this->dist[s];
    auto &parent = this->parent[s];

    dist[v] += diff;

    for (auto [to, cost] : G[v]) {
        if (parent[to] != v) continue;
        dfs_sub_tree(s, to, diff);
    }
}

void graph::dijkstra_remove_edge(edge e) {
    auto [v0, v1, cost] = e;

    auto v0_target = find(G[v0].begin(), G[v0].end(), std::make_pair(v1, cost));
    G[v0].erase(v0_target);
    // G[v0].shrink_to_fit();
    
    auto v1_target = find(G[v1].begin(), G[v1].end(), std::make_pair(v0, cost));
    G[v1].erase(v1_target);
    // G[v1].shrink_to_fit();

    remove_edge_que.emplace(e);

    for (auto& [s, dist] : dist) {
        auto &parent = this->parent[s];

        pque q;

        if (v0 == parent[v1]) {
            // v0 -> v1
            dfs(s, v1, e, q);
        }
        if (v1 == parent[v0]) {
            // v1 -> v0
            dfs(s, v0, e, q);
        }

        dijkstra_sub_func(s, q);
    }
}

long long graph::sum() {
    long long res = 0;

    for (auto& [s, dist] : dist) {
        res += std::accumulate(dist.begin(), dist.end(), 0ll);
    }

    return res;
}

void graph::backup() {
    for (const auto& [key, value] : dist) {
        dist_bk[key] = value; // 手動でコピーする
    }
    for (const auto& [key, value] : parent) {
        parent_bk[key] = value; // 手動でコピーする
    }
}

void graph::rollback() {
    for (const auto& [key, value] : dist_bk) {
        dist[key] = value; // 手動でコピーする
    }
    for (const auto& [key, value] : parent_bk) {
        parent[key] = value; // 手動でコピーする
    }

    while (add_edge_que.size()) {
        auto [from, to, cost] = add_edge_que.front();
        add_edge_que.pop();

        auto from_target = find(G[from].begin(), G[from].end(), std::make_pair(to, cost));
        G[from].erase(from_target);
        G[from].shrink_to_fit();
        
        auto to_target = find(G[to].begin(), G[to].end(), std::make_pair(from, cost));
        G[to].erase(to_target);
        G[to].shrink_to_fit();
    }

    while (remove_edge_que.size()) {
        auto [from, to, cost] = remove_edge_que.front();
        remove_edge_que.pop();

        add_edge(from, to, cost);
    }

    commit();
}

void graph::commit() {
    add_edge_que = std::queue<edge>();
    remove_edge_que = std::queue<edge>();
    dist_bk.clear();
    parent_bk.clear();
}