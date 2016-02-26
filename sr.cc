#include "generator/generator.h"
#include "ranklookup.h"
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <sstream>
#include <queue>
#include <vector>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

template <class RankLookupType>
class SRSat {
public:
    void read(std::string filename);
    void create(std::vector<std::vector<int> > instance);  // create from generated instance
    void shrink();
    void show();
    bool phase2();
    std::vector<int> find_rotation(std::vector<int>& fst, std::vector<int>& snd,
            std::vector<int>& last, int x);
    SRSat();
    ~SRSat();
//protected:
//    void set_rank(int i, int j, int rank);
//    int get_rank(int i, int j);

private:
    RankLookupType rank_lookup;

    int n;

    std::vector<std::vector<int> > pref; // Same as in Patrick's SR.java:
                // pref[i][k] = j <-> agent_i has agent_j as k^th choice

    std::vector<int> length; // length of agent's preference list

};

template<class RankLookupType>
SRSat<RankLookupType>::SRSat() {
}

template<class RankLookupType>
SRSat<RankLookupType>::~SRSat() {
}

/*
 * reads an instance from a file
 */
template<class RankLookupType>
void SRSat<RankLookupType>::read(std::string filename) {
    std::ifstream file;
    file.open(filename);

    std::string line;
    std::getline(file, line);
    std::istringstream is(line);

    is >> n;
    
    rank_lookup.initialise(n);

    for (int i=0; i<n; i++) {
        pref.push_back(std::vector<int>());
    }

    for (int i=0; i<n; i++) {
        std::getline(file, line);
        is.str(line);
        is.clear();

        int j;
        int k=0;
        while (is >> j) {
            rank_lookup.set_rank(i, j-1, k);
            pref[i].push_back(j - 1);
            ++k;
        }
        length.push_back(k);
        rank_lookup.reserve_space(i, k+1);
        rank_lookup.set_rank(i, i, k);
        pref[i].push_back(i);
    }

    file.close();
}

template<class RankLookupType>
void SRSat<RankLookupType>::create(std::vector<std::vector<int> > instance) {
    n = instance.size();

    rank_lookup.initialise(n);

    for (int i=0; i<n; i++) {
        pref.push_back(instance[i]);
    }

    for (int i=0; i<n; i++) {
        length.push_back(instance[i].size());
        rank_lookup.reserve_space(i, instance[i].size()+1);
        for (int j=0; j<instance[i].size(); j++) {
            int k = instance[i][j];
            rank_lookup.set_rank(i, k, j);
        }
        rank_lookup.set_rank(i, i, length[i]);
        pref[i].push_back(i);
    }
}

int digit_count(int x) {
    int dc = 0;
    while (x > 0) {
        ++dc;
        x /= 10;
    }
    return dc;
}

template<class RankLookupType>
void SRSat<RankLookupType>::show() {
    std::streamsize ss = std::cout.width();
    int w = digit_count(n);
    for (int i=0; i<n; i++) {
            std::cout.width(w);
            std::cout << (i+1) << " :";
        for (int j=0; j<length[i]; j++) {
            std::cout.width(j ? w+2 : w+1);
            std::cout << (pref[i][j]+1);
        }
        std::cout << std::endl;
    }
    std::cout.width(ss);

}

class Pair {
public:
    int a;
    int b;
    Pair(int a, int b): a(a), b(b) {}
};


/*
 * Reduce the problem instance size
 */
template<class RankLookupType>
void SRSat<RankLookupType>::shrink() {
    // For each agent i in q, we look at the person j that q likes most,
    // and remove everything worse than agent i in j's preference list
    std::queue<Pair> q;

    // For each of the n agents, we keep track of the index of the first
    // non-removed item on their preference list
    std::vector<int> first(n, 0);

    // Can agent i be matched with his jth preference?
    std::vector<std::vector<bool> > possible;
    for (int i=0; i<n; i++) {
        // TODO: Maybe the +1 isn't needed here?
        possible.push_back(std::vector<bool>(length[i]+1, true));
    }

    for (int i=0; i<n; i++) {
        if (length[i] > 0) {
            q.push(Pair(i, pref[i][0]));
        }
    }
    
    while (!q.empty()) {
        Pair p = q.front();
        q.pop();
        int j = p.b;
        int k = rank_lookup.get_rank(j, p.a, pref); // k is position of i in j's pref list
        for (int l=k+1; l<length[j]; l++) { // l is a position in j's pref list
            int m = pref[j][l];  // m is someone who ranks worse than i in
                                 // j's preference list
            if (possible[j][l]) {
                possible[j][l] = false;
                int rmj = rank_lookup.get_rank(m, j, pref);
                possible[m][rmj] = false;
                if (rmj == first[m]) {
                    do {
                        first[m]++;
                        // TODO: I think that this if statement could be moved
                        // below the do-while loop, or something like that?
                        if (pref[m][first[m]] != m) {
                            q.push(Pair(m, pref[m][first[m]]));
                        }
                    } while (!possible[m][first[m]]);
                }
            }
        }
        if (length[j] > k+1) length[j] = k+1;
    }

    rank_lookup.clear();
    for (int i=0; i<n; i++) {
        int j;   // other agent
        int k=0; // position in shrunk preference list
        for (int l=0; l<length[i]; l++) {
            j=pref[i][l];
            if (possible[i][l]) {
                rank_lookup.set_rank(i, j, k);
                pref[i][k] = j;
                ++k;
            }
        }
        rank_lookup.set_rank(i, i, k);
        pref[i][k] = i;
        length[i] = k;
    }
    
}


///////////////////////////////////////////
////////  Phase 2
///////////////////////////////////////////

template<class RankLookupType>
std::vector<int> SRSat<RankLookupType>::find_rotation(
        std::vector<int>& fst, std::vector<int>& snd, std::vector<int>& last, int x) {
    // Find a rotation, returning the x agent from each pair
    // x is the current x agent in the walk

    std::vector<int> walk;
    std::vector<int> first_pos_in_walk(n, -1);
    
    while (first_pos_in_walk[x] == -1) {
//        std::cout << x << " " << pref[x][fst[x]] << "!" << std::endl;
        first_pos_in_walk[x] = walk.size();
        walk.push_back(x);
        int y = pref[x][snd[x]]; // x's second-preference agent
        x = pref[y][last[y]];
    }
    // TODO: walk around graph

    std::vector<int> rotation(walk.begin() + first_pos_in_walk[x], walk.end());

    return rotation;
}

template<class RankLookupType>
bool SRSat<RankLookupType>::phase2() {
    // Can agent i be matched with his jth preference?
    std::vector<std::vector<bool> > possible;
    // How many elements of possible[i] are set to true?
    std::vector<int> num_possible(n);
    int num_at_least_2 = 0; // Number of agents with at least 2 possible agents on pref list
    int first_with_at_least_2 = n; // First agent with at least 2 possible agents on pref list
    for (int i=0; i<n; i++) {
        possible.push_back(std::vector<bool>(length[i], true));
        num_possible[i] = length[i];
        if (num_possible[i] >= 2) {
            num_at_least_2++;
            if (first_with_at_least_2 == n) first_with_at_least_2 = i;
        }
    }

    // fst[i] is the position in i's list of i's first preference.
    // snd[i] is the position in i's list of i's second preference.
    // These are -1 if no such element exists.
    std::vector<int> fst(n);
    std::vector<int> snd(n);
    std::vector<int> last(n);
    for (int i=0; i<n; i++) {
        // TODO: maybe we don't need to set -1? -- just check last[i]?
        if (length[i] == 0) fst[i] = -1; // Otherwise it's zero
        snd[i] = length[i]==1 ? -1 : 1;
        last[i] = length[i]-1;
        //std::cout << "    " << pref[i][length[i]] << std::endl;
    }
    
    while (first_with_at_least_2 < n) {
        // remove rotation
        std::vector<int> rotation_xs = find_rotation(fst, snd, last, first_with_at_least_2);
        std::vector<int> rotation_ynexts(rotation_xs.size());
        for (int i=0; i<rotation_xs.size(); i++) {
            int x = rotation_xs[i];
            rotation_ynexts[i] = pref[x][snd[x]]; // y of next pair in rotation
        }

        for (int i=0; i<rotation_xs.size(); i++) {
            int x = rotation_xs[i];
            int ynext = rotation_ynexts[i];
            //std::cout << x << " " << ynext << " " << std::endl;
            int rx = rank_lookup.get_rank(ynext, x, pref); // rank of x in ynext's pref list
            //std::cout << rx << "#" << std::endl;
            for (int j=rx+1; j<=last[ynext]; j++) {
                if (possible[ynext][j]) {
                    num_possible[ynext]--;
                    possible[ynext][j] = false;  // <-- TODO unnecessary, I think.
                    
                    // delete ynext from pref list of pref[ynext][j]
                    int k = pref[ynext][j];
                    int pos = rank_lookup.get_rank(k, ynext, pref);
                    num_possible[k]--;
                    if (num_possible[k] == 0) {
                        //std::cout << "-- UNSTABLE --" << std::endl;
                        return false;
                    }
                    possible[k][pos] = false;
                    if (pos==fst[k]) {
                        if (num_possible[k] >= 1) {
                            fst[k] = snd[k];
                            //std::cout << "---" << fst[k] << snd[k] << std::endl;
                        }
                        if (num_possible[k] >= 2) {
                            do {
                                ++snd[k];
                            } while (!possible[k][snd[k]]);
                        }
                    } else if (pos==snd[k] && num_possible[k] >= 2) {
                        // set new snd[k]
                        while (!possible[k][snd[k]]) ++snd[k];
                    }
                }
            }
            last[ynext] = rx;
        }
//        for (int i=0; i<n; i++) {
//            for (int j=0; j<length[i]; j++) {
//                std::cout.width(4);
//                std::cout << possible[i][j];
//            }
//            std::cout << std::endl;
//            std::cout << fst[i] << " " << snd[i] << " " << last[i] << "   " << num_possible[i] << std::endl;
//            std::cout << std::endl;
//        }


        while (num_possible[first_with_at_least_2] < 2)
            first_with_at_least_2++;

        //std::cout << "Press Enter to Continue";
        //std::cin.ignore();
    }

//    std::cout << "-- STABLE --" << std::endl;

//    for (int i=0; i<n; i++) {
//        if (num_possible[i])
//            std::cout << i << " " << pref[i][fst[i]] << std::endl;
//    }
    return true;
}

///////////////////////////////////////////
////////  End of Phase 2
///////////////////////////////////////////

template<class RankLookupType>
int normal_run(std::string filename, bool verbose) {
    clock_t start_time;
    clock_t shrink_start_time;
    clock_t solve_start_time;
    double read_time;
    double shrink_time;
    double solve_time;
    double total_time;

    start_time = clock();

    std::ios_base::sync_with_stdio(false);

    SRSat<RankLookupType> srSat;

    srSat.read(filename);
    read_time = double(clock() - start_time)/CLOCKS_PER_SEC;

    shrink_start_time = clock();
    srSat.shrink();
    shrink_time = double(clock() - shrink_start_time)/CLOCKS_PER_SEC;

    if (verbose) srSat.show();

    solve_start_time = clock();
    bool is_stable = srSat.phase2();
    solve_time = double(clock() - solve_start_time)/CLOCKS_PER_SEC;
   
    total_time = double(clock() - start_time)/CLOCKS_PER_SEC;

    std::cout << "Result:       " << (is_stable ? "STABLE" : "UNSTABLE") << std::endl;
    std::cout << "Read time:    " << read_time << std::endl;
    std::cout << "Shrink time:  " << shrink_time << std::endl;
    std::cout << "Phase 2 time: " << solve_time << std::endl;
    std::cout << "Total time:   " << total_time << std::endl;

	return 0;
}

template<class RankLookupType>
int random_run(double timeout, unsigned int n, double p, unsigned int seed) {

    std::ios_base::sync_with_stdio(false);

    if (p > 1) {
        std::cout << n << "\t" << p << "\t" << 0 << "\t" << 
                    -1 << "\t" << seed << std::endl;

        return 0;
    }

    int stable_count = 0;
    unsigned int num_instances = 0;
    clock_t start_time = clock();
    std::mt19937 rgen;
    rgen.seed(seed);
    while (true) {
        SRSat<RankLookupType> srSat;

        srSat.create(generate(n, p, rgen));

        srSat.shrink();

        num_instances++;

        bool is_stable = srSat.phase2();

        if (is_stable) stable_count++;

        if (double(clock() - start_time)/CLOCKS_PER_SEC > timeout) break;
    }

    double proportion_stable = (double)stable_count / (double)num_instances;
    std::cout << n << "\t" << p << "\t" << num_instances << "\t" << 
                proportion_stable << "\t" << seed << std::endl;

    return 0;
}

int main(int argc, char** argv) {
    double timeout = -1;
    unsigned int n = 100;
    double p;
    unsigned int seed;
    int type;
    bool verbose;
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "show help message")
            ("file,f", po::value<std::string>(), "read the specified file")
            ("random,r", "create a random instance")
            ("timeout", po::value<double>(&timeout)->default_value(5),
                    "the number of seconds after which to stop generating instances")
            ("n,n", po::value<unsigned int>(&n)->default_value(100), "number of agents")
            ("p,p", po::value<double>(&p)->default_value(0), "p")
            ("seed,s", po::value<unsigned int>(&seed)->default_value(1),
                    "seed for random generator")
            ("type", po::value<int>(&type)->default_value(1),
                    "1=array, 2=map, 3=linear scan")
            ("verbose,v", po::bool_switch(&verbose),
                    "Display verbose output (this is not fully implemented yet)")
            ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (type < 1 || type > 3) {
            std::cout << "Invalid type." << std::endl << desc << std::endl;
            return 1;
        }
        if (vm.count("help") || !(vm.count("file") || vm.count("random"))) {
            std::cout << desc << "\n";
            return 1;
        } else if (vm.count("file")) {
            std::string filename = vm["file"].as<std::string>();
            switch (type) {
                case 1: return normal_run<RankLookupArray>(filename, verbose);
                case 2: return normal_run<RankLookupMap>(filename, verbose);
                case 3: return normal_run<RankLookupLinearScan>(filename, verbose);
            }
        } else if (vm.count("random")) {
            switch (type) {
                case 1: return random_run<RankLookupArray>(timeout, n, p, seed);
                case 2: return random_run<RankLookupMap>(timeout, n, p, seed);
                case 3: return random_run<RankLookupLinearScan>(timeout, n, p, seed);

            }
        }
    } catch (po::error& e) {
        std::cout << e.what() << std::endl;
        return 1;
    }
}

