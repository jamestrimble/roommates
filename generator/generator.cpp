#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <cmath>
#include <sstream>
#include <utility>
#include <vector>

#include "generator.h"

inline int how_many_to_choose(double p,
                     unsigned int t, // maximum possible number to choose
                     std::mt19937_64& rgen) {
    std::binomial_distribution<int> binom_dist(t, p);
    return binom_dist(rgen);
}


/**
 * Generates a random subset of integers in the range [0, end),
 * with each value being chosen independently with probability p.
 * big_list is a vector of integers 0, 1, ..., end-1. At the end
 * of the method, the numbers [0, ..., end) might be in the wrong
 * order in big_list
 */
std::vector<int> random_subset(int end, double p,
                std::vector<int>& big_list, std::mt19937_64& rgen) {
    int how_many = how_many_to_choose(p, end, rgen);

    std::vector<int> result;
    if (how_many > 0) {

        for (int i=0; i<how_many; i++) {
            // Do a partial Knuth Shuffle
            std::uniform_int_distribution<int> gen(i, end-1);
            int r = gen(rgen);    
            int temp = big_list[i];
            big_list[i] = big_list[r];
            big_list[r] = temp;
        }

        std::copy(big_list.begin(), big_list.begin() + how_many, std::back_inserter(result));
    }
    return result;
}

std::vector<std::vector<int> > generate(unsigned int n, double p, std::mt19937_64& rgen) {
    std::vector<int> big_list;
    for (unsigned int i=0; i<n; i++) {
        big_list.push_back(i);
    }

    std::vector<std::vector<int> > pref_lists;

    // The first agent's preference list is initially empty
    pref_lists.push_back(std::vector<int>());
    
    // for each agent i (i>0), find which agents j with j<i
    // i ranks
    for (unsigned int i=1; i<n; i++) {
        pref_lists.push_back(random_subset(i, p, big_list, rgen));
        for (auto j : pref_lists[i]) {
            pref_lists[j].push_back(i);
        }
    }

    for (unsigned int i=0; i<n; i++) {
        std::vector<int>& v = pref_lists[i];
        std::shuffle(v.begin(), v.end(), rgen);
    }

    return pref_lists;
}

