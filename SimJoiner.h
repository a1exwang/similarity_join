#ifndef __EXP2_SIMJOINER_H__
#define __EXP2_SIMJOINER_H__

#include <vector>
#include <functional>
#include <vector>
#include <map>
#include <fstream>

template <typename IDType, typename SimType>
struct JoinResult {
    IDType id1;
    IDType id2;
    SimType s;
};

typedef JoinResult<unsigned, double> JaccardJoinResult;
typedef JoinResult<unsigned, unsigned> EDJoinResult;

const int SUCCESS = 0;
const int FAILURE = 1;

inline int editDist1(const char *str1 , const char *str2 , uint32_t m , uint32_t n);

class SimJoiner {
public:
    SimJoiner();
    ~SimJoiner();

    int joinJaccard(const char *filename1, const char *filename2, double threshold, std::vector<JaccardJoinResult> &result);
    int joinED(const char *filename1, const char *filename2, unsigned tau, std::vector<EDJoinResult> &result);

private:
     void createInvertedIndex(
        std::map<std::string, std::vector<int>> &idx,
        const std::vector<std::vector<std::string>> &lines
    );
    void createInvertedIndex(
        std::map<std::string, std::vector<int>> &idx,
        const std::vector<std::vector<std::string>> &lines,
        std::function<int (int len)> cb
    );

    void readFile(std::vector<std::vector<std::string>> &lines, std::ifstream &fs);
    static double jaccard(const std::vector<std::string> &s1,
                          const std::vector<std::string> &s2,
                          std::function<bool (const std::string &, const std::string &)>);
    static void sortWordsByIDF(std::vector<std::string> &words, const std::map<std::string, int> &dict);
    static std::function<bool (const std::string &, const std::string &)> getIDFComparator(
        const std::map<std::string, int> &dict
    );
private:
    // ED
    static int editDist(const std::string &s1, const std::string &s2);

    void readFile(std::vector<std::string> &lines, std::ifstream &fs);
    void createEDInvertedIndex(
        std::map<int, std::map<int, std::map<std::string, std::vector<int>>>> &idx,
        const std::vector<std::string> &lines,
        int tau,
        std::function<void (int l, int rid, int p, const std::string &segment)> cb);

    void createEDInvertedIndex(
        std::map<int, std::map<int, std::map<std::string, std::vector<int>>>> &idx,
        const std::vector<std::string> &lines,
        int tau);
};

#endif
