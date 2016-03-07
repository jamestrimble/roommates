#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <cmath>
#include <sstream>
#include <utility>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include "generator.h"

std::vector<std::vector<int> > generate_morph(unsigned int n, double p, std::mt19937_64& rgen) {
    // Using Patrick's recipe
    // TODO: throw an error if n is an odd number

    // TODO: Refactor to reduce duplication of shuffle

    boost::dynamic_bitset<> E1(n*n);  // As adjacency matrix. For j>i, E1[i*n + j] == 1 <-> edge (i,j) exists

    std::vector<std::pair<int, int> > clique_edges;
    for (int i=0; i<n-1; i++)
        for (int j=i+1; j<n; j++)
            clique_edges.push_back(std::pair<int, int>(i, j));

    int sz = n*n / 4;
    for (int i=0; i<sz; i++) {
        // Do a partial Knuth Shuffle, adding each element to the adjacency matrix rather than the array
        std::uniform_int_distribution<int> gen(i, clique_edges.size()-1);
        int r = gen(rgen);
        auto edge = clique_edges[r];
        E1[edge.first*n + edge.second] = true;
        clique_edges[r] = clique_edges[i];
    }

    boost::dynamic_bitset<> E2(n*n);  // As adjacency matrix. For j>i, E2[i*n + j] == 1 <-> edge (i,j) exists
    for (int i=0; i<n/2; i++)
        for (int j=n/2; j<n; j++)
            E2[i*n + j] = true;

    boost::dynamic_bitset<> Ex = E1 & E2;

    boost::dynamic_bitset<> E1minusEx = E1 - Ex;
    boost::dynamic_bitset<> E2minusEx = E2 - Ex;

    std::vector<std::pair<int, int> > E1minusEx_edges;
    std::vector<std::pair<int, int> > E2minusEx_edges;
    for (int i=0; i<n-1; i++) {
        for (int j=i+1; j<n; j++) {
            if (E1minusEx[i*n + j]) E1minusEx_edges.push_back(std::pair<int, int>(i, j));
            if (E2minusEx[i*n + j]) E2minusEx_edges.push_back(std::pair<int, int>(i, j));
        }
    }

    boost::dynamic_bitset<> chosen_edges = Ex;

    sz = p * E1minusEx_edges.size();
    for (int i=0; i<sz; i++) {
        // Do a partial Knuth Shuffle, adding each element to the adjacency matrix rather than the array
        std::uniform_int_distribution<int> gen(i, E1minusEx_edges.size()-1);
        int r = gen(rgen);
        auto edge = E1minusEx_edges[r];
        chosen_edges[edge.first*n + edge.second] = true;
        E1minusEx_edges[r] = E1minusEx_edges[i];
    }

    sz = (1-p) * E2minusEx_edges.size();
    for (int i=0; i<sz; i++) {
        // Do a partial Knuth Shuffle, adding each element to the adjacency matrix rather than the array
        std::uniform_int_distribution<int> gen(i, E2minusEx_edges.size()-1);
        int r = gen(rgen);
        auto edge = E2minusEx_edges[r];
        chosen_edges[edge.first*n + edge.second] = true;
        E2minusEx_edges[r] = E2minusEx_edges[i];
    }

    std::vector<std::vector<int> > pref_lists;
    for (int i=0; i<n; i++) pref_lists.push_back(std::vector<int>());

    for (int i=0; i<n-1; i++) {
        for (int j=i+1; j<n; j++) {
            if (chosen_edges[i*n + j]) {
                pref_lists[i].push_back(j);
                pref_lists[j].push_back(i);
            }
        }
    }

    for (unsigned int i=0; i<n; i++) {
        std::vector<int>& v = pref_lists[i];
        std::shuffle(v.begin(), v.end(), rgen);
    }

    return pref_lists;
}

