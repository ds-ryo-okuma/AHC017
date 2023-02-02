#include <iostream>
#include <vector>
#include <set>
#include <numeric>
#include <climits>

unsigned int randxor() {
    static unsigned int x = 123456789, y = 362436069, z = 521288629, w = 88675123;
    unsigned int t;
    t = (x ^ (x << 11)); x = y; y = z; z = w; return(w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
}

static double rand01() {
    return (randxor() + 0.5) * (1.0 / UINT_MAX);
}

clock_t base_ms = 0;
void tic() { base_ms = clock(); };
clock_t toc() { return clock() - base_ms; };

long long comp(const std::vector<long long> &A, const std::vector<long long> &B) {
    return accumulate(A.begin(), A.end(), 0ll) - accumulate(B.begin(), B.end(), 0ll);
}

void output (int M, const std::vector<std::set<int>> &S) {
    std::vector<int> R(M);
    for (int d = 0; d < (int)S.size(); ++d) {
        for (auto &i : S[d]) {
            R[i] = d;
        }
    }
    for (auto &r : R) std::cout << r + 1 << " ";
}