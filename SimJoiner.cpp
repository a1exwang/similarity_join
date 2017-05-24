#include "SimJoiner.h"
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <algorithm>


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
    if (i1 <= i2) {
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
  }

  return SUCCESS;
}

int SimJoiner::joinED(const char *filename1, const char *filename2, unsigned threshold, vector<EDJoinResult> &result) {
  result.clear();
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
