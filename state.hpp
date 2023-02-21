#pragma once

#include <vector>

#include "util.hpp"
#include "graph.hpp"

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