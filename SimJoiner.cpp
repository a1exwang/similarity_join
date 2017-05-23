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

  for (auto l : idx1) {
    cout << l.first << ": ";
    for (auto i : l.second) {
      cout << i << ' ';
    }
    cout << endl;
  }

  // I. prefix filter
  auto prefixFn = [threshold](int x) -> int {
      return int((1 - threshold) * x + 1);
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

  // III. suffix positional filter

  // IV. final filter

  for (auto p : pairs) {
    double j = jaccard(lines1[p.first], lines2[p.second]);
    if (j >= threshold) {
      result.push_back({p.first, p.second, j});

      cout << "file1[" << p.first << "]( ";
      for (auto w : lines1[p.first]) {
        cout << w << ' ';
      }
      cout << "), " << "file2[" << p.second << "]( ";
      for (auto w : lines2[p.second]) {
        cout << w << ' ';
      }
      cout << "), " << j << endl;
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

double SimJoiner::jaccard(const std::vector<std::string> &s1, const std::vector<std::string> &s2) {
  vector<string> inters;
  vector<string> uni;
  set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(inters));
  set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(uni));
  return (double)inters.size() / uni.size();
}
