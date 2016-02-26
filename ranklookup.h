#include <vector>
#include <unordered_map>

using std::vector;

class RankLookupLinearScan {
public:
    inline void set_rank(int i, int j, int rank) {};
    inline int get_rank(int i, int j, vector<vector<int> >& pref) {
        for (vector<vector<int> >::size_type k=0; k<pref[i].size(); k++)
            if (pref[i][k]==j)
                return k;
        return -1; // should never be reached
    }
    void initialise(int n) {};
    void reserve_space(int i, int size) {};
    void clear(int i) {};
    ~RankLookupLinearScan() {};
};


class RankLookupMap {
public:
    inline void set_rank(int i, int j, int rank) { this->rank[i][j] = rank; }
    inline int get_rank(int i, int j, vector<vector<int> >& pref) {
        return rank[i][j];
    }
    void initialise(int n) {
        for (int i=0; i<n; i++)
            rank.push_back(std::unordered_map<int, int>());
    }
    void reserve_space(int i, int size) { rank[i].reserve(size); }
    void clear(int i) {
         rank[i].clear();
    }
    ~RankLookupMap() {}
private:
    vector<std::unordered_map<int, int> > rank;
};


class RankLookupArray {
public:
    inline void set_rank(int i, int j, int rank) { this->rank[i][j] = rank; }
    inline int get_rank(int i, int j, vector<vector<int> >& pref) {
        return rank[i][j];
    }
    void initialise(int n) {
        rank = new int*[n];
        int* rankBlock = new int[n*n];
        for (int i=0; i<n; i++)
            rank[i] = rankBlock + i*n;
    }
    void reserve_space(int i, int size) {}
    void clear(int i) {}
    ~RankLookupArray() {
        if (rank) {
            delete[] rank[0];
            delete[] rank;
        }
    }
private:
    // THIS IS ASSIGNED AS A SINGLE BLOCK OF MEMORY
    //vector<vector<int> > rank; // Same as in Patrick's SR.java:
    int** rank; // Same as in Patrick's SR.java:
                // rank[i][j] = k <-> agent_i ranks agent_j as k^th choice
    
};

