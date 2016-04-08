#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <cmath>
#include <sstream>
#include <utility>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "Generator.h"

void shuffle_pref_lists(std::vector<std::vector<int> >& pref_lists, int n, std::mt19937_64& rgen) {
    for (int i=0; i<n; i++) {
        std::vector<int>& v = pref_lists[i];
        std::shuffle(v.begin(), v.end(), rgen);
    }
}

///////////////////////////////////////////////////////////////

// https://networkx.github.io/documentation/latest/_modules/networkx/generators/random_graphs.html#fast_gnp_random_graph
std::vector<std::vector<int> > GeneratorEdgeGeneration::generate() {
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    
    std::vector<std::vector<int> > pref_lists(n, std::vector<int>());

    int w = -1;
    double lp = std::log(1.0 - p);

    int v = 1;
    while (v < n) {
        double lr = std::log(1.0 - dis(rgen));
        if (p > 0.999999999)
            ++w;
        else
            w = w + 1 + (int)(lr/lp);
        while (w >= v && v < n) {
            w = w - v;
            ++v;
        }
        if (v < n) {
            pref_lists[v].push_back(w);
            pref_lists[w].push_back(v);
        }
    }

    shuffle_pref_lists(pref_lists, n, rgen);
    return pref_lists;
}

///////////////////////////////////////////////////////////////

std::vector<std::vector<int> > GeneratorEdgeGenerationSimple::generate() {
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    
    std::vector<std::vector<int> > pref_lists(n, std::vector<int>());

    for (int i=0; i<n-1; i++) {
        for (int j=i+1; j<n; j++) {
            if (dis(rgen) < p) {
                pref_lists[i].push_back(j);
                pref_lists[j].push_back(i);
            }
        }
    }

    shuffle_pref_lists(pref_lists, n, rgen);
    return pref_lists;
}

///////////////////////////////////////////////////////////////

inline int how_many_to_choose(double p,
                     int t, // maximum possible number to choose
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

std::vector<std::vector<int> > GeneratorEdgeGenerationBinom::generate() {
    std::vector<int> big_list;
    for (int i=0; i<n; i++) {
        big_list.push_back(i);
    }

    std::vector<std::vector<int> > pref_lists;

    // The first agent's preference list is initially empty
    pref_lists.push_back(std::vector<int>());
    
    // for each agent i (i>0), find which agents j with j<i
    // i ranks
    for (int i=1; i<n; i++) {
        pref_lists.push_back(random_subset(i, p, big_list, rgen));
        for (auto j : pref_lists[i]) {
            pref_lists[j].push_back(i);
        }
    }

    shuffle_pref_lists(pref_lists, n, rgen);
    return pref_lists;
}

///////////////////////////////////////////////////////////////

std::vector<std::vector<int> > GeneratorEdgeGenerationComplement::generate() {
    std::vector<std::vector<int> > complement(comp_gen.generate());

    std::vector<std::vector<int> > pref_lists(n, std::vector<int>());

    boost::dynamic_bitset<> comp_adj_mat(n*n);
    for (int i=0; i<n; i++)
        for (auto j : complement[i])
            comp_adj_mat[i*n + j] = 1;

    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            if (!comp_adj_mat[i*n + j] && i != j)
                pref_lists[i].push_back(j);

    shuffle_pref_lists(pref_lists, n, rgen);
    return pref_lists;
}

///////////////////////////////////////////////////////////////

void shuffle_to_adjmat(std::vector<std::pair<int, int> >& edges, boost::dynamic_bitset<>& adj_mat,
        int m, int n, std::mt19937_64& rgen) {
    // Select m edges from _edges_ by shuffling, and add them to adj_mat.
    for (int i=0; i<m; i++) {
        // Do a partial Knuth Shuffle, adding each element to the adjacency matrix rather than the array
        std::uniform_int_distribution<int> gen(i, edges.size()-1);
        int r = gen(rgen);
        auto edge = edges[r];
        adj_mat[edge.first*n + edge.second] = true;
        std::swap(edges[r], edges[i]);
    }
}

std::vector<std::vector<int> > GeneratorSMMorph::generate() {
    // Using Patrick's recipe
    
    if (n % 1 != 0) throw GenError("n must be even for this generator");

    boost::dynamic_bitset<> E1(n*n);  // As adjacency matrix. For j>i, E1[i*n + j] == 1 <-> edge (i,j) exists

    std::vector<std::pair<int, int> > clique_edges;
    for (int i=0; i<n-1; i++)
        for (int j=i+1; j<n; j++)
            clique_edges.push_back(std::pair<int, int>(i, j));

    shuffle_to_adjmat(clique_edges, E1, (n/2)*(n/2), n, rgen);

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

    unsigned int m1 = std::round(p * E1minusEx_edges.size());
    if (m1 > E1minusEx_edges.size()) m1 = E1minusEx_edges.size();
    shuffle_to_adjmat(E1minusEx_edges, chosen_edges, m1, n, rgen);

    unsigned int m2 = std::round((1-p) * E2minusEx_edges.size());
    if (m2 > E2minusEx_edges.size()) m2 = E2minusEx_edges.size();
    shuffle_to_adjmat(E2minusEx_edges, chosen_edges, m2, n, rgen);

    std::vector<std::vector<int> > pref_lists(n, std::vector<int>());

    for (int i=0; i<n-1; i++) {
        for (int j=i+1; j<n; j++) {
            if (chosen_edges[i*n + j]) {
                pref_lists[i].push_back(j);
                pref_lists[j].push_back(i);
            }
        }
    }

    shuffle_pref_lists(pref_lists, n, rgen);
    return pref_lists;
}

std::vector<std::vector<int> > GeneratorSMMorphTypeA::generate() {
    if (n % 1 != 0) throw GenError("n must be even for this generator");

    int m = (n/2) * (n/2);

    boost::dynamic_bitset<> E1(n*n);  // As adjacency matrix. For j>i, E1[i*n + j] == 1 <-> edge (i,j) exists
    for (int i=0; i<n-1; i++)
        for (int j=i+1; j<n; j++)
            E1[i*n + j] = true;

    std::vector<std::pair<int, int> > biclique_edges;
    for (int i=0; i<n/2; i++)
        for (int j=n/2; j<n; j++)
            biclique_edges.push_back(std::pair<int, int>(i, j));

    boost::dynamic_bitset<> E(n*n);

    int m2 = std::round((1-p) * m);
    if (m2 > m) m2 = m;
    shuffle_to_adjmat(biclique_edges, E, m2, n, rgen);

    boost::dynamic_bitset<> E1minusE = E1 - E;

    std::vector<std::pair<int, int> > E1minusE_edges;
    for (int i=0; i<n-1; i++)
        for (int j=i+1; j<n; j++)
            if (E1minusE[i*n + j]) E1minusE_edges.push_back(std::pair<int, int>(i, j));

    int m1 = std::round(p * m);
    if (m1 > m) m1 = m;
    shuffle_to_adjmat(E1minusE_edges, E, m1, n, rgen);

    std::vector<std::vector<int> > pref_lists(n, std::vector<int>());

    for (int i=0; i<n-1; i++) {
        for (int j=i+1; j<n; j++) {
            if (E[i*n + j]) {
                pref_lists[i].push_back(j);
                pref_lists[j].push_back(i);
            }
        }
    }

    shuffle_pref_lists(pref_lists, n, rgen);
    return pref_lists;
}

///////////////////////////////////////////////////////////////

std::vector<std::vector<int> > GeneratorEdgeSelection::generate() {
    int m = std::round(p * n * (n-1) / 2.0); // number of edges
    if (m > n * (n-1) / 2) m = n * (n-1) / 2; // TODO: is this necessary?

    std::vector<std::vector<int> > pref_lists(n, std::vector<int>());

    for (int i=0; i<m; i++) {
        // Do a partial Knuth Shuffle
        std::uniform_int_distribution<int> gen(i, possible_edges.size() - 1);
        int r = gen(rgen);
        std::swap(possible_edges[i], possible_edges[r]);
        int v = possible_edges[i].first;
        int w = possible_edges[i].second;
        pref_lists[v].push_back(w);
        pref_lists[w].push_back(v);
    }

    shuffle_pref_lists(pref_lists, n, rgen);
    return pref_lists;
}

///////////////////////////////////////////////////////////////

std::vector<std::vector<int> > GeneratorCompleteGraph::generate() {
    if (p != 1) throw GenError("p must equal 1 for this generator");

    std::vector<std::vector<int> > pref_lists(n, std::vector<int>());

    for (int i=0; i<n-1; i++) {
        for (int j=i+1; j<n; j++) {
            pref_lists[i].push_back(j);
            pref_lists[j].push_back(i);
        }
    }

    shuffle_pref_lists(pref_lists, n, rgen);
    return pref_lists;
}
