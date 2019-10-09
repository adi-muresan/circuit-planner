#include "disjoint_sets.h"

utils::disj_sets::disj_sets(int size)
  : size(size), sets(size), rank(size, 0) {
  // each set is independent i.e. its own representative
  for(int i = 0; i < size; ++ i) {
    sets[i] = i;
  }
}

int utils::disj_sets::get_representative(int id) {
  // path compression
  int ancestor_id = sets[id];
  if(ancestor_id != id) {
    ancestor_id = get_representative(ancestor_id);
    sets[id] = ancestor_id;
  }
  return ancestor_id;
}

void utils::disj_sets::merge(int id1, int id2) {
  link(get_representative(id1), get_representative(id2));
}

void utils::disj_sets::link(int id1, int id2) {
  // union by rank
  if(rank[id1] > rank[id2]) {
    sets[id2] = id1;
  } else {
    sets[id1] = id2;
    if(rank[id1] == rank[id2]) {
      rank[id2] ++;
    }
  }
}
