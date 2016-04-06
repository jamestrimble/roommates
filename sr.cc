#include "generator/Generator.h"
#include "ranklookup.h"
#include "Solver.h"
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <sstream>
#include <queue>
#include <vector>

#include <boost/program_options.hpp>

using namespace Minisat;
namespace po = boost::program_options;

struct SolvingError: public std::runtime_error {
    SolvingError(std::string const& message)
        : std::runtime_error(message)
    {}
};

template <class RankLookupType>
class SRSat {
public:
    void read(std::string filename);
    void create(std::vector<std::vector<int> > instance);  // create from generated instance
    void phase1();
    void show();
    bool phase2();
    void build_sat();
    int solve_sat(bool display_sols);  // Find all solutions, and return count
    bool is_sat();
    int matching_size;
    SRSat();
    ~SRSat();

private:
    RankLookupType rank_lookup;
    vec<Lit> lits;  // SAT literals
    int n;   // number of agents
    std::vector<std::vector<int> > pref; // Same as in Patrick's SR.java:
                // pref[i][k] = j <-> agent_i has agent_j as k^th choice
    std::vector<int> length; // length of agent's preference list
    Solver s;   // MiniSat solver

    // --- MiniSat variables ---
    // assigned[i][j] <-> agent_i gets her j^th choice
    std::vector<std::vector<Lit> > assigned;
    // MiniSat variables are just integers. Keep track of the next one to be created.
    int n_vars;
    // Create a new variable
    inline Lit createLit();
    // Create an implies constraint
    inline void implies(Lit a, Lit b) {s.addClause(~a, b);}
    void displayMatching();
    std::vector<int> find_rotation(std::vector<int>& fst, std::vector<int>& snd,
            std::vector<int>& last, int x);
    void calc_matching_size();
};

template<class RankLookupType>
SRSat<RankLookupType>::SRSat(): n_vars(0) {
    s.ccmin_mode = 0;
    s.phase_saving = 0;
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

    for (int i=0; i<n; i++)
        pref.push_back(std::vector<int>());

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

    pref = instance;   // TODO: avoid copying the instance

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
    while (x > 0) { ++dc; x /= 10; }
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

template<class RankLookupType>
void SRSat<RankLookupType>::phase1() {
    // This implementation is partly based on Mertens' paper on random instances.

    // q is the list of free agents.
    // For each agent i in q, we look at the person j that i likes most,
    // and remove everything worse than agent i in j's preference list
    std::queue<int> q;

    // For each of the n agents, we keep track of the index of the first
    // non-removed item on their preference list
    std::vector<int> first(n, 0);
    
    // If prop_from[i] == j, then j is semiengaged to i
    std::vector<int> prop_from(n, -1);

    for (int i=0; i<n; i++)
        if (length[i] > 0)
            q.push(i);
    
    process_queue:
    while (!q.empty()) {
        int i = q.front();
        q.pop();

        int j = pref[i][first[i]];
        int k = rank_lookup.get_rank(j, i, pref); // k is position of i in j's pref list
        while (k >= length[j]) {
            if (i == j) goto process_queue;
            first[i]++;
            j = pref[i][first[i]];
            k = rank_lookup.get_rank(j, i, pref);
        }

        if (prop_from[j] != -1)
            q.push(prop_from[j]);   // old proposer is freed; add him to queue

        prop_from[j] = i;

        if (length[j] > k+1) length[j] = k+1;
    }

    std::vector<std::vector<int> > shrunk_pref;

    for (int i=0; i<n; i++) {
        shrunk_pref.push_back(std::vector<int>());
        for (int l=first[i]; l<length[i]; l++) {
            int j=pref[i][l];
            if (rank_lookup.get_rank(j, i, pref) < length[j])
                shrunk_pref[i].push_back(j);
        }
        shrunk_pref[i].push_back(i);
    }

    rank_lookup.clear();
    for (int i=0; i<n; i++) {
        for (std::vector<int>::size_type k=0; k<shrunk_pref[i].size(); k++)
            rank_lookup.set_rank(i, shrunk_pref[i][k], k);

        length[i] = shrunk_pref[i].size() - 1;
    }
    pref = std::move(shrunk_pref);
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
    // first_pos_in_walk[i] is the first position in the walk at which i
    // appears, or -1 if i does not appear in the walk.
    std::vector<int> first_pos_in_walk(n, -1);
    
    while (first_pos_in_walk[x] == -1) {
        first_pos_in_walk[x] = walk.size();
        walk.push_back(x);
        int y = pref[x][snd[x]]; // x's second-preference agent
        x = pref[y][last[y]];
    }

    // return the rotation (the walk without the tail)
    return std::vector<int>(walk.begin() + first_pos_in_walk[x], walk.end());
}

template<class RankLookupType>
bool SRSat<RankLookupType>::phase2() {
    // possible[i] == true <-> agent i can be matched with his jth preference
    std::vector<std::vector<bool> > possible;
    for (int i=0; i<n; i++) possible.push_back(std::vector<bool>(length[i], true));

    int first_with_at_least_2 = n; // First agent with at least 2 possible agents on pref list

    // fst[i] is the position in i's list of i's first preference.
    // snd[i] is the position in i's list of i's second preference.
    // last[i] is the position in i's list of i's last preference.
    std::vector<int> fst(n);
    std::vector<int> snd(n, 1);
    std::vector<int> last(n);
    for (int i=0; i<n; i++) {
        if (length[i] == 0) fst[i] = -1; // Otherwise it's zero
        last[i] = length[i]-1;
        if (length[i] >= 2 && first_with_at_least_2 == n)
            first_with_at_least_2 = i;
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
            int rx = rank_lookup.get_rank(ynext, x, pref); // rank of x in ynext's pref list
            for (int j=rx+1; j<=last[ynext]; j++) {
                if (possible[ynext][j]) {
                    //possible[ynext][j] = false;  // <-- TODO unnecessary, I think.
                    
                    // delete ynext from pref list of pref[ynext][j]
                    int k = pref[ynext][j];
                    int pos = rank_lookup.get_rank(k, ynext, pref);
                    if (fst[k] == last[k]) {
                        return false;
                    } else if (pos==last[k]) {
                        // TODO: can this ever happen?
                        std::cout << "eek!" << std::endl;
                    }

                    possible[k][pos] = false;
                    if (pos==fst[k]) {
                        fst[k] = snd[k];
                        if (fst[k] != last[k])
                            do { ++snd[k]; } while (!possible[k][snd[k]]);
                    } else if (pos==snd[k]) {
                        // set new snd[k]
                        while (!possible[k][snd[k]]) ++snd[k];
                    }
                }
            }
            last[ynext] = rx;
        }

        // Update our 'pointer' to the first agent with at least 2 remaining preferences
        while (fst[first_with_at_least_2] == last[first_with_at_least_2])
            first_with_at_least_2++;
    }

    return true;
}

///////////////////////////////////////////
////////  End of Phase 2
///////////////////////////////////////////


///////////////////////////////////////////
////////  SAT stuff
///////////////////////////////////////////

/*
 * Creates a literal representing a new MiniSat variable
 */
template<class RankLookupType>
Lit SRSat<RankLookupType>::createLit() {
    Lit lit = mkLit(s.newVar());
    lits.push(lit);
    n_vars++;
    return lit;
}

/*
 * Builds the SAT instance
 */
template<class RankLookupType>
void SRSat<RankLookupType>::build_sat() {
    for (int i=0; i<n; i++) {
        assigned.push_back(std::vector<Lit>(length[i]+1, mkLit(0)));
    }
    std::vector<std::vector<Lit> > assignedBelow(n, std::vector<Lit>());
    std::vector<std::vector<Lit> > assignedAbove(n, std::vector<Lit>());

    //Lit* assignedBlock = new Lit[n*n];
    //for (int i=0; i<n; i++) {
    //    assigned[i] = assignedBlock + i*n;
    //}

    // For each agent A, create a variable which is true iff A
    // is unassigned
    for (int i=0; i<n; i++) {
        assigned[i][length[i]] = createLit();       
    }

    for (int i=0; i<n; i++) {
        // For each agent A, create a variable for each position in A's preference
        // list which is true iff A gets better than her j^th choice. Similarly,
        // create variables which are true iff A gets worse than her j^th choice.
        for (int j=0; j<=length[i]; j++) {
            assignedBelow[i].push_back(createLit());
            assignedAbove[i].push_back(createLit());
        }
    }

    // Create variables, each of which is true iff agent is is matched with
    // agent k. Do this only if i<k. Below, we will re-use these variables for the
    // case where i>k.
    for (int i=0; i<n; i++) {
        for (int j=0; j<length[i]; j++) {
            int k = pref[i][j];
            if (i<k) assigned[i][j] = createLit();
        }
    }

    // Use the same variable for agent A being matched with agent B as for agent B
    // being matched with agent A.
    for (int i=0; i<n; i++) {
        for (int j=0; j<length[i]; j++) {
            int k = pref[i][j];
            if (i>k) assigned[i][j] = assigned[k][rank_lookup.get_rank(k, i, pref)];
        }
    }

    for (int i=0; i<n; i++) {
        // Constraint: For each agent A, either A is unmatched (second term), or
        // is matched to someone else (second term)
        s.addClause(assignedBelow[i][length[i]], assigned[i][length[i]]);
    }

    for (int i=0; i<n; i++) {
        // Each agent can't do worse than being assigned to herself
        s.addClause(~assignedAbove[i][length[i]]);
        // Each agent can't do better than her first preference
        s.addClause(~assignedBelow[i][0]);
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<length[i]; j++) {
            int k = pref[i][j];
            // Equivalent to the following constraint in SR.java:
            // implies(gt(agent[i],rank[i][k]),lt(agent[k],rank[k][i]))
            implies(assignedAbove[i][rank_lookup.get_rank(i, k, pref)],
                    assignedBelow[k][rank_lookup.get_rank(k, i, pref)]);

            // It isn't necessary to have as many variables and
            // constraints as below, but this system seems to work
            // quite well. I haven't yet tested a version similar
            // to the model of Prosser and Gent for stable marriage

            // assignedAbove[i][j+1] implies assignedAbove[i][j]
            implies(assignedAbove[i][j+1], assignedAbove[i][j]);
            // assigned[i][j+1] implies assignedAbove[i][j]
            implies(assigned[i][j+1], assignedAbove[i][j]);
            // Agent i can't be assigned to her j^th preference and also
            // assigned to someone worse than her j^th preference
            implies(assignedAbove[i][j], ~assigned[i][j]);
            // If agent i is assigned to worse than her j^th preference,
            // then either she is assigned to her (j+1)^th preference or
            // she is assigned to worse than her (j+1)^th preference
            s.addClause(~assignedAbove[i][j],
                        assigned[i][j+1],
                        assignedAbove[i][j+1]);

            implies(assignedBelow[i][j], assignedBelow[i][j+1]);
            implies(assigned[i][j], assignedBelow[i][j+1]);
            implies(assignedBelow[i][j+1], ~assigned[i][j+1]);
            s.addClause(~assignedBelow[i][j+1],
                        assigned[i][j],
                        assignedBelow[i][j]);
        }
    }
}

template<class RankLookupType>
void SRSat<RankLookupType>::displayMatching() {
    for (int i=0; i<n; i++) {
        for (int j=0; j<length[i]; j++) {
            int k = pref[i][j];
            if (i<k && s.modelValue(assigned[i][j]) == l_True) {
                std::cout << "(" << i+1 << "," << k+1 << ") ";
            }
        }
    }
    std::cout << std::endl;
}

template<class RankLookupType>
void SRSat<RankLookupType>::calc_matching_size() {
    matching_size = 0;
    for (int i=0; i<n; i++) {
        for (int j=0; j<length[i]; j++) {
            int k = pref[i][j];
            if (i<k && s.modelValue(assigned[i][j]) == l_True) matching_size += 2;
        }
    }
}

template<class RankLookupType>
int SRSat<RankLookupType>::solve_sat(bool display_sols) {
    // Return value: number of stable solutions
    int n_solutions = 0;

    vec<Lit> newClause;
    while (s.solve()) {
        if (display_sols) displayMatching();
        n_solutions++;
        if (n_solutions == 1) calc_matching_size();

        // Create a new clause to rule out found solution
        newClause.clear();
        for (int i=0; i<lits.size(); i++) {
            if (s.modelValue(lits[i])!=l_True)
                newClause.push(lits[i]);
            else
                newClause.push(~lits[i]);
        }

        s.addClause(newClause);
    }

    return n_solutions;
}

template<class RankLookupType>
bool SRSat<RankLookupType>::is_sat() {
    bool stable_sol_exists = s.solve();
    if (stable_sol_exists) calc_matching_size();
    return stable_sol_exists;
}


///////////////////////////////////////////
////////  End of SAT stuff
///////////////////////////////////////////


template<class RankLookupType>
int normal_run(std::string filename, bool verbose, bool all_sols) {
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
    srSat.phase1();
    shrink_time = double(clock() - shrink_start_time)/CLOCKS_PER_SEC;

    if (verbose) srSat.show();

    solve_start_time = clock();
    bool is_stable = srSat.phase2();

    if (all_sols && is_stable) {
        srSat.build_sat();
        int n_sols = srSat.solve_sat(true);
        if (n_sols)
            std::cout << n_sols << "\tstable solutions found" << std::endl;
    }

    solve_time = double(clock() - solve_start_time)/CLOCKS_PER_SEC;
   
    total_time = double(clock() - start_time)/CLOCKS_PER_SEC;

    std::cout << "Result:       " << (is_stable ? "STABLE" : "UNSTABLE") << std::endl;
    std::cout << "Read time:    " << read_time << std::endl;
    std::cout << "Phase 1 time:  " << shrink_time << std::endl;
    std::cout << "Phase 2 time: " << solve_time << std::endl;
    std::cout << "Total time:   " << total_time << std::endl;

	return 0;
}

template<class RankLookupType, class GeneratorType>
int do_random_run(double timeout, int max_iter, int n, double p, bool all_sols,
               std::mt19937_64& rgen, bool use_phase_1, int seed, bool record_sol_sizes) {

    std::ios_base::sync_with_stdio(false);

    if (p > 1) {
        // TODO: throw an exception
        return -1;
    }

    auto gen = GeneratorType(n, p, rgen);

    int stable_count = 0;
    int num_instances = 0;
    clock_t start_time = clock();

    // A map with stable solution count as key and number of instances
    // with this number of stable sols as value
    // This map is only used if all_sols is true
    std::map<int, int> sol_count_map;
    if (all_sols) sol_count_map[0] = 0;

    // A map with number of agents in a stable solution as key and number of
    // instances with this many agents in a stable solution as value.
    // This map is only used if record_sol_sizes is true
    std::map<int, int> sol_size_map;

    while (true) {
        SRSat<RankLookupType> srSat;
        srSat.create(gen.generate());
        if (use_phase_1) srSat.phase1();
        num_instances++;
        srSat.build_sat();
        bool is_stable;
        if (all_sols)  {
            int n_sols = srSat.solve_sat(false);
            if (n_sols > 0) stable_count++;
            if (sol_count_map.count(n_sols))
                ++sol_count_map[n_sols];
            else
                sol_count_map[n_sols] = 1;
            if (use_phase_1 && (n_sols==0) == srSat.phase2()) throw SolvingError("Solvers disagree!");
            is_stable = n_sols > 0;
        } else {
            is_stable = srSat.is_sat();
            if (is_stable) stable_count++;
            // Double-check correctness
            if (use_phase_1 && is_stable != srSat.phase2()) throw SolvingError("Solvers disagree!");
        }
        if (is_stable && record_sol_sizes) {
            int sol_size = srSat.matching_size;
            if (sol_size_map.count(sol_size))
                ++sol_size_map[sol_size];
            else
                sol_size_map[sol_size] = 1;
        }
        if ((max_iter!=-1 && num_instances==max_iter) || double(clock() - start_time)/CLOCKS_PER_SEC > timeout) break;
    }

    if (all_sols)
        for (auto it : sol_count_map)
            std::cout << "sol_count\t" << it.first << "\t" << it.second << std::endl;
    if (record_sol_sizes)
        for (auto it : sol_size_map)
            std::cout << "sol_size\t" << it.first << "\t" << it.second << std::endl;
    
    double proportion_stable = (double)stable_count / (double)num_instances;
    std::cout << n << "\t" << p << "\t" << num_instances << "\t" << 
                proportion_stable << "\t" << seed << std::endl;

    return 0;
}

template<class RankLookupType>
int random_run(double timeout, int max_iter, int n, double p, bool all_sols,
               std::mt19937_64& rgen, bool use_phase_1, int seed, int gen_type, bool record_sol_sizes) {
    if (gen_type == 8) {
        if      (p == 1)  gen_type = 7;
        else if (p > 0.7) gen_type = 4;
        else              gen_type = 2;
    }
    switch(gen_type) {
    case 1:
        return do_random_run<RankLookupType, GeneratorEdgeGenerationSimple>(timeout, max_iter, n, p, all_sols,
                rgen, use_phase_1, seed, record_sol_sizes);
    case 2:
        return do_random_run<RankLookupType, GeneratorEdgeGeneration>(timeout, max_iter, n, p, all_sols,
                rgen, use_phase_1, seed, record_sol_sizes);
    case 3:
        return do_random_run<RankLookupType, GeneratorEdgeGenerationBinom>(timeout, max_iter, n, p, all_sols,
                rgen, use_phase_1, seed, record_sol_sizes);
    case 4:
        return do_random_run<RankLookupType, GeneratorEdgeGenerationComplement>(timeout, max_iter, n, p, all_sols,
                rgen, use_phase_1, seed, record_sol_sizes);
    case 5:
        return do_random_run<RankLookupType, GeneratorEdgeSelection>(timeout, max_iter, n, p, all_sols,
                rgen, use_phase_1, seed, record_sol_sizes);
    case 6:
        return do_random_run<RankLookupType, GeneratorSMMorph>(timeout, max_iter, n, p, all_sols,
                rgen, use_phase_1, seed, record_sol_sizes);
    case 7:
        return do_random_run<RankLookupType, GeneratorCompleteGraph>(timeout, max_iter, n, p, all_sols,
                rgen, use_phase_1, seed, record_sol_sizes);
    }
    return 1;
}

int main(int argc, char** argv) {
    double timeout = -1;
    int n = 100;
    double p, np;
    int seed;
    int type;
    int max_iter;
    bool verbose;
    bool all_sols;  // find all solutions, or just check whether stable?
    bool record_sol_sizes;  // record size of stable solutions? (random runs only)
    int gen_type;
    bool no_phase_1;
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "show help message")
            ("file,f", po::value<std::string>(), "read the specified file")
            ("random,r", "create a random instance")
            ("maxiter", po::value<int>(&max_iter)->default_value(-1), "maximum number of random runs")
            ("timeout", po::value<double>(&timeout)->default_value(5),
                    "the number of seconds after which to stop generating instances")
            ("n,n", po::value<int>(&n)->default_value(100), "number of agents")
            ("p,p", po::value<double>(&p)->default_value(0.5), "p (use --p or --np, not both)")
            ("np", po::value<double>(&np), "n*p (use --p or --np, not both)")
            ("seed,s", po::value<int>(&seed)->default_value(1),
                    "seed for random generator")
            ("type", po::value<int>(&type)->default_value(1),
                    "1=array, 2=map, 3=linear scan, 4=auto")
            ("verbose,v", po::bool_switch(&verbose),
                    "Display verbose output (this is not fully implemented yet)")
            ("all,a", po::bool_switch(&all_sols),
                    "Find all solutions? (Default: just check whether a stable solution exists)")
            ("record-sizes", po::bool_switch(&record_sol_sizes),
                    "Record sizes of stable solutions in random runs? (Default: off)")
            ("gen-type", po::value<int>(&gen_type)->default_value(2),
                     "generator type: 1=simple edge gen, 2=fast edge gen, "
                     "3=edge gen (using binomial dist),"
                     "4=edge gen (using complement) 5=edge selection, 6=SR/SM morph"
                     "7=complete, 8=auto")
            ("no-phase-1", po::bool_switch(&no_phase_1),
                    "Don't carry out phase 1? (Uses SAT solver only) (Ignored if input is from file)")
            ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (type < 1 || type > 4) {
            std::cout << "Invalid type." << std::endl << desc << std::endl;
            return 1;
        }

        if (gen_type < 1 || gen_type > 8) {
            std::cout << "Invalid generator type." << std::endl << desc << std::endl;
            return 1;
        }

        if (vm.count("np")) {
            p = np / n;
        }
        if (p < 0 || p > 1) {
            std::cout << "Invalid value for p." << std::endl << desc << std::endl;
            return 1;
        }
        
        if (type == 4) {
            if (n > 5000 || (n > 1000 && n*p < 200) || (n > 500 && n*p < 100))
                type = 3;
            else
                type = 1;
        }

        if (vm.count("help") || !(vm.count("file") || vm.count("random"))) {
            std::cout << desc << "\n";
            return 1;
        } else if (vm.count("file")) {
            std::string filename = vm["file"].as<std::string>();
            switch (type) {
                case 1: return normal_run<RankLookupArray>(filename, verbose, all_sols);
                case 2: return normal_run<RankLookupMap>(filename, verbose, all_sols);
                case 3: return normal_run<RankLookupLinearScan>(filename, verbose, all_sols);
            }
        } else if (vm.count("random")) {
            std::mt19937_64 rgen;
            rgen.seed(seed);

            // TODO: avoid passing seed
            switch (type) {
                case 1: 
                    return random_run<RankLookupArray>(
                            timeout, max_iter, n, p, all_sols, rgen, !no_phase_1, seed, gen_type, record_sol_sizes);
                case 2:
                    return random_run<RankLookupMap>(
                            timeout, max_iter, n, p, all_sols, rgen, !no_phase_1, seed, gen_type, record_sol_sizes);
                case 3:
                    return random_run<RankLookupLinearScan>(
                            timeout, max_iter, n, p, all_sols, rgen, !no_phase_1, seed, gen_type, record_sol_sizes);
            }
        }
    } catch (po::error& e) {
        std::cout << e.what() << std::endl;
        return 1;
    } catch (SolvingError& e) {
        std::cout << e.what() << std::endl;
        return 1;
    }
}

