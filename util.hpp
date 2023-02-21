#pragma once

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