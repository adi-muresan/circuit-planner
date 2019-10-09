#ifndef DISJOINT_SETS_H
#define DISJOINT_SETS_H

#include <vector>

namespace utils {

/* Implementation of the disjoint sets data structure from CLR */
class disj_sets {
  int size;
  std::vector<int> sets;
  std::vector<int> rank;

  void link(int id1, int id2);

public:
  disj_sets(int size);

  int get_representative(int id);

  void merge(int id1, int id2);
};

}

#endif // DISJOINT_SETS_H
