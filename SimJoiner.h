#ifndef __EXP2_SIMJOINER_H__
#define __EXP2_SIMJOINER_H__

#include <vector>
#include <functional>
#include <vector>
#include <map>
#include <fstream>
#include <set>

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

class SimJoiner {
public:
    SimJoiner();
    ~SimJoiner();

    int joinJaccard(const char *filename1, const char *filename2, double threshold, std::vector<JaccardJoinResult> &result);
    int joinED(const char *filename1, const char *filename2, unsigned threshold, std::vector<EDJoinResult> &result);

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
    static double jaccard(const std::vector<std::string> &s1, const std::vector<std::string> &s2);
};

#endif
