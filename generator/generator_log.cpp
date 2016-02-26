#include <algorithm>
#include <iostream>
#include <random>
#include <cmath>
#include <sstream>
#include <vector>

#include "generator.h"

inline int rand_skip(double log_1_minus_p,
                     std::uniform_real_distribution<double>& dis,
                     std::mt19937_64& rgen) {
    double r = dis(rgen);
    if (r==0) r = 0.000000000000000000001;
    return (int) (std::log(r) / log_1_minus_p);
}

std::vector<std::vector<int> > generate(unsigned int n, double p, std::mt19937_64& rgen) {
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    std::vector<std::vector<int> > pref_lists(n, std::vector<int>());

    double log_1_minus_p = std::log(1-p);

    for (int i=0; i<n-1; i++) {
        int j = i;
        while(true) {
            j += rand_skip(log_1_minus_p, dis, rgen) + 1;
            if (j >= n) break;
            pref_lists[i].push_back(j);
            pref_lists[j].push_back(i);
        };
    }

    for (int i=0; i<n; i++) {
        std::vector<int>& v = pref_lists[i];
        std::shuffle(v.begin(), v.end(), rgen);
    }

    return pref_lists;
}

