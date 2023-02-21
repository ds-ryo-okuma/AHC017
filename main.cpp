#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <numeric>
#include <cmath>

#include "util.hpp"
#include "state.hpp"
#include "graph.hpp"
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
