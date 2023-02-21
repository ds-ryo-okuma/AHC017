#line 1 "main.cpp"
#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <numeric>
#include <cmath>

#line 2 "util.hpp"

#line 7 "util.hpp"
#include <climits>

unsigned int randxor() {
    static unsigned int x = 123456789, y = 362436069, z = 521288629, w = 88675123;
    unsigned int t;
    t = (x ^ (x << 11)); x = y; y = z; z = w; return(w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
}

static double rand01() {
    return (randxor() + 0.5) * (1.0 / UINT_MAX);
}

std::pair<int, int> rand2(int N) {
    int i = randxor() % N, j = randxor() % (N - 1);
    if (i <= j) j++;
    return {i, j};
}

clock_t base_ms = 0;
void tic() { base_ms = clock(); };
clock_t toc() { return clock() - base_ms; };

template <class T>
T comp(const std::vector<T> &A, const std::vector<T> &B) {
    return accumulate(A.begin(), A.end(), 0.0) - accumulate(B.begin(), B.end(), 0.0);
}

template <class T>
void output (int M, const std::vector<T> &S) {
    std::vector<int> R(M);
    for (int d = 0; d < (int)S.size(); ++d) {
        for (auto &i : S[d]) {
            R[i] = d;
        }
    }
    for (auto &r : R) std::cout << r + 1 << " ";
}
#line 2 "state.hpp"

#line 4 "state.hpp"

#line 2 "graph.hpp"

#line 5 "graph.hpp"
#include <unordered_map>
#line 7 "graph.hpp"
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

    for (auto& [s, dist] : this->dist) {
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

    for (auto& [s, dist] : this->dist) {
        res += std::accumulate(dist.begin(), dist.end(), 0ll);
    }

    return res;
}

void graph::backup() {
    for (const auto& [key, value] : this->dist) {
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
#line 7 "state.hpp"

const int START_NODES_NUM = 1;

class state {
public:
    // 頂点数
    int N;

    // 辺数
    int M;

    // 期間
    int D;

    // 1 日当たりの最大工事数
    int K;

    // 工事計画 R[i] := i 番目の辺を何日目に工事するか
    std::vector<int> R;

    // 工事予定数 num[i] := i 日目の工事数
    std::vector<int> num;

    // 辺集合
    std::vector<edge> edges;

    // グラフ集合 graphs[i] := i 日目に工事する辺を除いたグラフ
    std::vector<graph> graphs;

    // 全ての辺を工事しないときのグラフ
    graph all_path;

    // 始点の集合
    std::vector<int> start_nodes;

    // スコア score[i] := i 日目の不満度の合計
    std::vector<long long> score;

    state(int N, int M, int D, int K);
    void add_edge(int id, int from, int to, int cost);
    void initialize();
    long long calc_score();
    void neighbor();
    void output();
};

/**
 * @brief Construct a new state::state object
 * 
 * @param N 頂点数
 * @param M 辺数
 * @param D 期間
 * @param K 1 日当たりの最大工事数
 */
state::state(int N, int M, int D, int K) : N(N), M(M), D(D), K(K), R(M), num(D), edges(M), graphs(D, graph(N, M)), all_path(N, M), score(D) {
    // 始点を START_NODES_NUM 個選ぶ
    for (int i = 0; i < START_NODES_NUM; ++i) {
        // start_nodes.emplace_back(randxor() % N);
        start_nodes.emplace_back(0);
    }
}

/**
 * @brief 辺を追加
 * 
 * @param id 
 * @param from 
 * @param to 
 * @param cost 
 */
void state::add_edge(int id, int from, int to, int cost) {
    edges[id] = edge(from, to, cost);
    all_path.add_edge(from, to, cost);
}

/**
 * @brief 初期化処理
 * 
 */
void state::initialize() {
    // 各始点からの最短経路を計算
    for (int s : start_nodes) {
        all_path.dijkstra(s);
    }

    // i 番目の辺を i % D 日目に工事する
    for (int i = 0; i < M; ++i) {
        R[i] = i % D; ++num[R[i]];
    }

    // 工事計画を反映したグラフを生成
    for (int i = 0; i < M; ++i) {
        for (int d = 0; d < D; ++d) {
            if (R[i] == d) continue;
            auto [from, to, cost] = edges[i];
            graphs[d].add_edge(from, to, cost);
        }
    }

    // スコアを計算
    for (int d = 0; d < D; ++d) {
        for (int s : start_nodes) {
            graphs[d].dijkstra(s);

            for (int g = 0; g < N; ++g) {
                score[d] += graphs[d].dist[s][g] - all_path.dist[s][g];
            }
        }
    }
}

void state::neighbor() {
    // 70% の確立で 2 辺を入れ替える
    if (randxor() % 100 < 70) {
        auto [e0, e1] = rand2(M);
        while (R[e0] == R[e1]) {
            auto [temp0, temp1] = rand2(M);
            e0 = temp0, e1 = temp1;
        }

        int day0 = R[e0], day1 = R[e1];

        graphs[day0].backup();
        graphs[day0].dijkstra_add_edge(edges[e0]);
        graphs[day0].dijkstra_remove_edge(edges[e1]);

        graphs[day1].backup();
        graphs[day1].dijkstra_add_edge(edges[e1]);
        graphs[day1].dijkstra_remove_edge(edges[e0]);

        long long day0_sum = 0, day1_sum = 0;

        for (int s : start_nodes) {
            for (int g = 0; g < N; ++g) {
                day0_sum += graphs[day0].dist[s][g] - all_path.dist[s][g];
                day1_sum += graphs[day1].dist[s][g] - all_path.dist[s][g];
            }
        }

        if (day0_sum + day1_sum < score[day0] + score[day1]) {
            graphs[day0].commit();
            graphs[day1].commit();

            score[day0] = day0_sum;
            score[day1] = day1_sum;

            std::swap(R[e0], R[e1]);
        } else {
            graphs[day0].rollback();
            graphs[day1].rollback();
        }
    }

    // 30% の確率で 1 辺を別の日に変更
    else {
        int e = randxor() % M;
        while (num[R[e]] == 0) {
            e = randxor() % M;
        }

        int day0 = randxor() % D;
        while (R[e] == day0 || num[day0] == K) {
            day0 = randxor() % D;
        }

        int day1 = R[e];

        graphs[day0].backup();
        graphs[day0].dijkstra_remove_edge(edges[e]);

        graphs[day1].backup();
        graphs[day1].dijkstra_add_edge(edges[e]);

        long long day0_sum = 0, day1_sum = 0;

        for (int s : start_nodes) {
            for (int g = 0; g < N; ++g) {
                day0_sum += graphs[day0].dist[s][g] - all_path.dist[s][g];
                day1_sum += graphs[day1].dist[s][g] - all_path.dist[s][g];
            }
        }

        if (day0_sum + day1_sum < score[day0] + score[day1]) {
            graphs[day0].commit();
            graphs[day1].commit();

            score[day0] = day0_sum;
            score[day1] = day1_sum;

            --num[day1];
            ++num[day0];
            R[e] = day1;
        } else {
            graphs[day0].rollback();
            graphs[day1].rollback();
        }
    }
}

void state::output() {
    for (int i = 0; i < M; ++i) {
        printf("%d%c", R[i] + 1, i < M - 1 ? ' ' : '\n');
    }
}
#line 11 "main.cpp"
using namespace std;


const double CALC_TIME = 5.8;

int main() {
    int N, M, D, K;
    cin >> N >> M >> D >> K;

    // 現在の状態
    state now(N, M, D, K);

    for (int i = 0; i < M; ++i) {
        int u, v, w;
        cin >> u >> v >> w;
        now.add_edge(i, --u, --v, w);
    }

    now.initialize();

    int loop_count = 0;

    tic();
    while (toc() < CALC_TIME * CLOCKS_PER_SEC) {
        now.neighbor();

        ++loop_count;
    }

    now.output();

    cerr << "INITIAL_COUNT = " << 0 << endl;
    cerr << "LOOP_COUNT = " << loop_count << endl;

    return 0;
}
