#include <vector>
#include <set>
using namespace std;

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
