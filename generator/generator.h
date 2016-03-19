#include <random>
#include <stdexcept>
#include <vector>

struct GenError: public std::runtime_error {
    GenError(std::string const& message)
        : std::runtime_error(message)
    {}
};

class Generator {
protected:
    unsigned int n;
    double p;
    std::mt19937_64& rgen;
public:
    Generator(unsigned int n, double p, std::mt19937_64& rgen) : n(n), p(p), rgen(rgen) {}
    virtual ~Generator() {}
    virtual std::vector<std::vector<int> > generate() = 0;
};

class GeneratorEdgeGeneration : public Generator {
public:
    GeneratorEdgeGeneration(unsigned int n, double p, std::mt19937_64& rgen)
        : Generator(n, p, rgen) {}
    std::vector<std::vector<int> > generate();
};

class GeneratorEdgeGenerationBinom: public Generator {
public:
    GeneratorEdgeGenerationBinom(unsigned int n, double p, std::mt19937_64& rgen)
        : Generator(n, p, rgen) {}
    std::vector<std::vector<int> > generate();
};

class GeneratorEdgeSelection: public Generator {
public:
    GeneratorEdgeSelection(unsigned int n, double p, std::mt19937_64& rgen)
        : Generator(n, p, rgen)
    {
        for (int i=0; i<n-1; i++)
            for (int j=i+1; j<n; j++)
                possible_edges.push_back(std::pair<int, int>(i, j));
    }
    std::vector<std::vector<int> > generate();
private:
    std::vector<std::pair<int, int> > possible_edges;
};

class GeneratorSMMorph : public Generator {
public:
    GeneratorSMMorph(unsigned int n, double p, std::mt19937_64& rgen)
        : Generator(n, p, rgen) {}
    std::vector<std::vector<int> > generate();
};
