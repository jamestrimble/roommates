#include <iostream>
#include <sstream>
#include <random>
#include <vector>

#include "Generator.h"

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cout << argv[0] << " n p gen_type seed" << std::endl;
        std::cout << "options for type: 1=simple edge gen, 2=fast edge gen, "
                     "3=edge gen (using binomial dist), "
                     "4=edge gen (using complement) 5=edge selection, 6=SR/SM morph type A, 7=SR/SM morph type B, "
                     "8=complete"
                     << std::endl;
    } else {
        int n;
        double p;
        int gen_type;
        int seed;

        std::istringstream iss(argv[1]);
        iss >> n;

        iss.str(argv[2]); iss.clear();
        iss >> p;

        iss.str(argv[3]); iss.clear();
        iss >> gen_type;

        iss.str(argv[4]); iss.clear();
        iss >> seed;

        std::mt19937_64 rgen;
        rgen.seed(seed);

        Generator * gen;
        switch(gen_type) {
            case 1: gen = new GeneratorEdgeGenerationSimple(n, p, rgen); break;
            case 2: gen = new GeneratorEdgeGeneration(n, p, rgen); break;
            case 3: gen = new GeneratorEdgeGenerationBinom(n, p, rgen); break;
            case 4: gen = new GeneratorEdgeGenerationComplement(n, p, rgen); break;
            case 5: gen = new GeneratorEdgeSelection(n, p, rgen); break;
            case 6: gen = new GeneratorSMMorphTypeA(n, p, rgen); break;
            case 7: gen = new GeneratorSMMorph(n, p, rgen); break;
            case 8: gen = new GeneratorCompleteGraph(n, p, rgen); break;
        }

        std::vector<std::vector<int> > pref_lists = gen->generate();

        std::cout << n << std::endl;
        for (int i=0; i<n; i++) {
            std::vector<int>& v = pref_lists[i];
            for (auto it=v.begin(); it!=v.end(); ++it)
                std::cout << (*it)+1 << " ";
            std::cout << std::endl;
        }

        delete gen;
    }
}


