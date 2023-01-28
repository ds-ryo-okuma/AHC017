#include <iostream>
#include <vector>
using namespace std;

int main() {
    int N, M, D, K;
    cin >> N >> M >> D >> K;

    vector<int> U(M), V(M), W(M);
    for (int i = 0; i < M; ++i) {
        cin >> U[i] >> V[i] >> W[i];
    }

    vector<int> X(N), Y(N);
    for (int i = 0; i < N; ++i) {
        cin >> X[i] >> Y[i];
    }

    vector<int> R(M);

    int day = 1, count = 0;
    for (auto &r : R) {
        if (count == K) {
            day++;
            count = 0;
        }
        r = day;
        count++;
    }

    for (auto &r : R) cout << r << " ";

    return 0;
}
