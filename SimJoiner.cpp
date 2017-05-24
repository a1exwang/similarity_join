#include "SimJoiner.h"
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <cassert>
#define MY_DEBUG


using namespace std;

SimJoiner::SimJoiner() {
}

SimJoiner::~SimJoiner() {
}

int SimJoiner::joinJaccard(const char *filename1,
                           const char *filename2,
                           double threshold,
                           vector<JaccardJoinResult> &result) {
  result.clear();
  // TODO
  // 1. threshold >= or >
  // 2. a a b ==~ a b c ?

  // Create inverted index
  ifstream f1(filename1);
  ifstream f2(filename2);
  if (!f1.is_open()) {
    cout << "bad";
  }
  vector<vector<string>> lines1;
  map<string, vector<int>> idx1;
  readFile(lines1, f1);
  createInvertedIndex(idx1, lines1);

  vector<vector<string>> lines2;
  map<string, vector<int>> idx2;
  readFile(lines2, f2);
  createInvertedIndex(idx2, lines2);

  map<string, int> idx0;
  for (auto &p : idx1) {
    if (idx0.find(p.first) == idx0.end()) {
      idx0[p.first] = (int)p.second.size();
    }
    else {
      idx0[p.first] += (int)p.second.size();
    }
  }
  for (auto &p : idx2) {
    if (idx0.find(p.first) == idx0.end()) {
      idx0[p.first] = (int)p.second.size();
    }
    else {
      idx0[p.first] += (int)p.second.size();
    }
  }

  for (auto &l : lines1)
    sortWordsByIDF(l, idx0);
  for (auto &l : lines2)
    sortWordsByIDF(l, idx0);

  // I. prefix filter
  auto prefixFn = [threshold](int x) -> int {
      return min({int((1 - threshold) * x + 1), x});
  };
  std::vector<std::pair<unsigned, unsigned>> pairs;
  map<string, vector<int>> idx1_prefix;
  map<string, vector<int>> idx2_prefix;

  //  1. create index for prefixes
  createInvertedIndex(idx1_prefix, lines1, prefixFn);
  createInvertedIndex(idx2_prefix, lines2, prefixFn);

  //  2. traverse indexes that has >= 2 lines, mark them valid
  for (auto idx : idx1_prefix) {
    auto &word = idx.first;
    auto &rids1 = idx.second;
    if (idx2_prefix.find(word) != idx2_prefix.end()) {
      auto &rids2 = idx2_prefix[word];

      for (auto rid1 : rids1) {
        for (auto rid2 : rids2) {
          pairs.push_back({rid1, rid2});
        }
      }
    }
  }

  // II. positional prefix filter
  std::vector<std::pair<unsigned, unsigned>> pairs2;
  for (auto p : pairs) {
    auto &l1 = lines1[p.first];
    auto &l2 = lines2[p.second];
    auto it11 = l1.begin();
    auto it12 = it11 + prefixFn((int)l1.size());

    auto it21 = l2.begin();
    auto it22 = it21 + prefixFn((int)l2.size());

    vector<string> intersection_set;
    set_intersection(it11, it12, it21, it22, back_inserter(intersection_set), getIDFComparator(idx0));

    vector<string> union_set;
    set_union(it11, it12, it21, it22, back_inserter(union_set), getIDFComparator(idx0));
    int ubound = (int)intersection_set.size() +
        std::min({(int)l1.size() - prefixFn((int)l1.size()), (int)l2.size() - prefixFn((int)l2.size())});
    int lbound = (int)union_set.size() +
                 std::max({(int)l1.size() - prefixFn((int)l1.size()), (int)l2.size() - prefixFn((int)l2.size())});
    if ((double)ubound / lbound >= threshold) {
      pairs2.push_back(p);
    }
  }

  // III. suffix positional filter

  // IV. final filter, and sort in the mean time
  map<pair<unsigned, unsigned>, double> result_set;
  for (auto p : pairs2) {
    double j = jaccard(lines1[p.first], lines2[p.second], getIDFComparator(idx0));
    // cout << p.first << ", " << p.second << ", " << j << endl;
    if (j >= threshold)
      result_set[{p.first, p.second}] = j;
  }

  // V. output
  for (auto r_item : result_set) {
    auto i1 = r_item.first.first;
    auto i2 = r_item.first.second;
    auto j = r_item.second;
    result.push_back({i1, i2, j});

//      cout << "file1[" << i1 << "]( ";
//      for (auto w : lines1[i1]) {
//        cout << w << ' ';
//      }
//      cout << "), " << "file2[" << i2 << "]( ";
//      for (auto w : lines2[i2]) {
//        cout << w << ' ';
//      }
//      cout << "), " << j << endl;
  }

  return SUCCESS;
}

int SimJoiner::joinED(const char *filename1, const char *filename2, unsigned tau, vector<EDJoinResult> &result) {
  result.clear();

  ifstream fs1(filename1);
  vector<string> lines1;
  ifstream fs2(filename2);
  vector<string> lines2;
  map<int, map<int, map<string, vector<int>>>> idx1;
  map<int, map<int, map<string, vector<int>>>> idx2;
  readFile(lines1, fs1);
  readFile(lines2, fs2);

  // 1. Create segment index for file2
  createEDInvertedIndex(idx2, lines2, tau);
//  for (auto a : idx2) {
//    for (auto b : a.second) {
//      for (auto c : b.second) {
//        for (auto i : c.second) {
//          cout << a.first << " " << b.first << " " << c.first << " " << i << endl;
//        }
//      }
//    }
//  }

  // 2. For each line in file1
  //  3. Segment the line, scan for the segments, and get candidates
  set<pair<int, int>> resultMap;
  for (int rid1 = 0; rid1 < lines1.size(); ++rid1) {
    auto &line1 = lines1[rid1];
    int l1 = (int) line1.size();
#ifdef MY_DEBUG
    cout << "rid1: " << rid1 << ' ' << line1 << endl;
#endif
    for (int deltaTau = -tau; deltaTau <= (int)tau; ++deltaTau) {
      int l2 = l1 + deltaTau;
      if (idx2.find(l2) != idx2.end()) {
#ifdef MY_DEBUG
        cout << "\tnearly_same_len: " << l2 << endl;
#endif
        int i = 1;
        for (auto pair2 : idx2[l2]) {
          // pos2 is pi
          int pos2 = pair2.first + 1;
          auto &segs = pair2.second;
          int delta = l1 - l2;
          int lowerBound = max({pos2 - (i - 1), pos2 + delta - ((int) tau + 1 - i)});
          int upperBound = min({pos2 + (i - 1), pos2 + delta + ((int) tau + 1 - i)});
#ifdef MY_DEBUG
          cout << "\t\tpos2 = " << pos2 << endl;
          cout << "\t\t[" << lowerBound << ", " << upperBound << "]" << endl;
#endif
          for (int startPos1 = lowerBound; startPos1 <= upperBound; ++startPos1) {
            auto segLen = segs.begin()->first.size();
            if (startPos1-1 + segLen <= line1.size()) {
              auto subs = line1.substr((unsigned)startPos1-1, segLen);
#ifdef MY_DEBUG
              cout << "\t\t\ttrying (" << startPos1 << ", " << subs << ")" << endl;
#endif
              if (segs.find(subs) != segs.end()) {
                const string &seg = subs;
                const vector<int> &ridOfSegs2 = segs[subs];
                for (auto rid2 : ridOfSegs2)
                  resultMap.insert({rid1, rid2});
              }
            }
          }
          i++;
        }
      }
    }
  }

#ifdef MY_DEBUG
  cout << "resultMap.size() = " << resultMap.size() << endl;
#endif

  // 4. Verification
  for (auto p : resultMap) {
    auto l1 = p.first;
    auto l2 = p.second;
    auto ed = editDist(lines1[l1], lines2[l2]);
    if (ed <= tau) {
      result.push_back({(unsigned)l1, (unsigned)l2, (unsigned)ed});
#ifdef MY_DEBUG
      cout << l1 << ", " << l2 << "   " << ed << endl;
#endif
    }
  }
  return SUCCESS;
}

void SimJoiner::createInvertedIndex(
    std::map<std::string, std::vector<int>> &idx,
    const std::vector<std::vector<std::string>> &lines,
    std::function<int(int len)> cb) {
  for (auto rid = 0; rid < (int)lines.size(); rid++) {
    auto &line = lines[rid];
    int words_needed = cb((int)lines.size());
    for (auto wid = 0; wid < words_needed; wid++) {
      auto &word = line[wid];
      if (idx.find(word) == idx.end()) {
        idx[word] = {rid};
      }
      else {
        idx[word].push_back(rid);
      }
    }
  }
}

void SimJoiner::createInvertedIndex(
    std::map<std::string, std::vector<int>> &idx,
    const std::vector<std::vector<std::string>> &lines) {
  auto cb = [](int x) -> int { return x; };
  createInvertedIndex(idx, lines, cb);
}

void split(const std::string &s, std::vector<std::string> &ret) {
  char delim = ' ';
  size_t last = 0;
  size_t index = s.find_first_of(delim,last);
  while (index!=std::string::npos) {
    ret.push_back(s.substr(last,index - last));
    last = index + 1;
    index = s.find_first_of(delim,last);
  }
  if (index-last>0) {
    ret.push_back(s.substr(last, index-last));
  }
}

void SimJoiner::readFile(vector<vector<string>> &lines, ifstream &fs) {
  string line;
  while (std::getline(fs, line)) {
    vector<string> words;
    split(line, words);
    lines.push_back(move(words));
  }
}

double SimJoiner::jaccard(const std::vector<std::string> &s1,
                          const std::vector<std::string> &s2,
                          std::function<bool (const std::string &, const std::string &)> cmp) {
  vector<string> inters;
  vector<string> uni;
  set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(inters), cmp);
  set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(uni), cmp);
  return (double)inters.size() / uni.size();
}

void SimJoiner::sortWordsByIDF(std::vector<std::string> &words, const std::map<std::string, int> &dict) {
  sort(words.begin(), words.end(), getIDFComparator(dict));
}

std::function<bool(const std::string &, const std::string &)>
SimJoiner::getIDFComparator(const std::map<std::string, int> &dict) {
  return [&dict](const string &lhs, const string &rhs) -> bool {
      if (dict.find(lhs)->second < dict.find(rhs)->second) {
        return true;
      }
      else if (dict.find(lhs)->second == dict.find(rhs)->second) {
        return lhs < rhs;
      }
      else {
        return false;
      }
  };
}

int SimJoiner::editDist(const std::string &s1, const string &s2) {
  int tmp[s1.size()+1][s2.size()+1];
  memset(tmp, 0, (s1.size()+1)*(s2.size()+1) * sizeof(int));

  int i = 0;
  int j = 0;
  while (true) {
    if (i == 0 && j == 0) {
      tmp[i][j] = 0;
    }
    else if (i == 0) {
      tmp[i][j] = tmp[i][j-1] + 1;
    }
    else if (j == 0) {
      tmp[i][j] = tmp[i-1][j] + 1;
    }
    else {
      int a = tmp[i][j - 1] + 1;
      int b = tmp[i-1][j] + 1;
      int c = tmp[i-1][j-1] + (s1[i - 1] == s2[j - 1] ? 0 : 1);
      tmp[i][j] = min({a, b, c});
    }

    if (i == (int)s1.size() && j == (int)s2.size()) {
      break;
    }

    if (j == 0 || i == (int)s1.size()) {
      if (i + j + 1 <= (int)s2.size()) {
        j = i + j + 1;
        i = 0;
      }
      else {
        i = i + j + 1 - (int)s2.size();
        j = (int)s2.size();
      }
    }
    else {
      i++;
      j--;
    }
  }
  auto result = tmp[s1.size()][s2.size()];
  return result;
}

void SimJoiner::readFile(std::vector<std::string> &lines, std::ifstream &fs) {
  string line;
  while (std::getline(fs, line)) {
    lines.push_back(line);
  }
}

void SimJoiner::createEDInvertedIndex(
    map<int, map<int, map<string, vector<int>>>> &idx,
    const std::vector<std::string> &lines,
    int tau) {


  for (int rid = 0; rid < lines.size(); rid++) {
    int l = (int)lines[rid].size();
    const string &line = lines[rid];
    if (idx.find(l) == idx.end()) {
      idx[l] = map<int, map<string , vector<int>>>();
    }
    map<int, map<string, vector<int>>> s2;
    map<string, vector<int>> s1;
    vector<int> segments;
    int segmentCount2 = l % (tau + 1);
    int segmentLen2 = (int)std::ceil((double)l / (tau + 1));
    int segmentCount1 = tau + 1 - segmentCount2;
    int segmentLen1 = (l - segmentCount2 * segmentLen2) / segmentCount1;
    assert(segmentLen1 == (int)std::floor((double) l / (tau + 1)));

    for (int j = 0; j < segmentCount1; ++j) {
      int p = j * segmentLen1;
      auto seg = line.substr((unsigned)p, (unsigned)segmentLen1);

      if (idx[l].find(p) == idx[l].end()) {
        idx[l][p] = map<string, vector<int>>();
      }
      if (idx[l][p].find(seg) == idx[l][p].end()) {
        idx[l][p][seg] = vector<int>();
      }
      idx[l][p][seg].push_back(rid);
    }
    for (int j = 0; j < segmentCount2; ++j) {
      int p = j * segmentLen2 + segmentCount1 * segmentLen1;
      auto seg = line.substr((unsigned)p, (unsigned)segmentLen2);

      if (idx[l].find(p) == idx[l].end()) {
        idx[l][p] = map<string, vector<int>>();
      }
      if (idx[l][p].find(seg) == idx[l][p].end()) {
        idx[l][p][seg] = vector<int>();
      }
      idx[l][p][seg].push_back(rid);
    }
  }
}

