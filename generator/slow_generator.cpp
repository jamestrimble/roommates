#include <algorithm>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

void generate(unsigned int n, double p, unsigned int seed) {
    std::cout << n << std::endl;

    std::mt19937 rgen;
    rgen.seed(seed);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    std::vector<int> v;
    std::vector<std::vector<bool> > A;
    for (int i=0; i<n; i++)
        A.push_back(std::vector<bool>(n, false));

    for (int i=0; i<n-1; i++)
        for (int j=i+1; j<n; j++)
            if (p >= dis(rgen))
                A[i][j] = A[j][i] = true;

    for (int i=0; i<n; i++) {
        v.clear();
        for (int j=0; j<n; j++) if (A[i][j]) v.push_back(j);
        std::shuffle(v.begin(), v.end(), rgen);
        for (auto it=v.begin(); it!=v.end(); ++it)
            std::cout << (*it)+1 << " ";
        std::cout << std::endl;
    }
}

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

        generate(n, p, seed);
    }
}

