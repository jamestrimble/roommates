#include <iostream>
#include <sstream>
#include <random>
#include <vector>

#include "generator.h"

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cout << argv[0] << " n p seed" << std::endl;
    } else {
        unsigned int n;
        double p;
        unsigned int seed;

        std::istringstream iss(argv[1]);
        iss >> n;

        iss.str(argv[2]); iss.clear();
        iss >> p;

        iss.str(argv[3]); iss.clear();
        iss >> seed;

        std::cout << n << std::endl;
        
        std::mt19937_64 rgen;
        rgen.seed(seed);
        std::vector<std::vector<int> > pref_lists = generate_morph(n, p, rgen);
        for (int i=0; i<n; i++) {
            std::vector<int>& v = pref_lists[i];
            for (auto it=v.begin(); it!=v.end(); ++it)
                std::cout << (*it)+1 << " ";
            std::cout << std::endl;
        }
    }
}


