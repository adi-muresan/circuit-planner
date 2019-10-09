#include "scoring.h"
#include "propagation.h"
#include "utils/disjoint_sets.h"

#include <cmath>
#include <limits>
#include <algorithm>

using namespace std;
using namespace propagation;

double scoring::compute_poly_distance(const std::vector<int> &target, const std::vector<int> &candidate) {
  if(candidate.empty()) {
    // infinite distance for nonexistent candidate
    return numeric_limits<double>::max();
  }

  double distance = 0;

  int best_match_id = 0;
  for(int p : target) {
    // find closest power in the candidate
    int pow_dist = abs(p - candidate[best_match_id]);

    while(best_match_id + 1 < candidate.size()
          && pow_dist > abs(p - candidate[best_match_id + 1])) {
      ++ best_match_id;
      pow_dist = abs(p - candidate[best_match_id]);
    }

    // TODO: there's a hidden hyperparameter here to represent the tradeoff / conversion
    distance += pow_dist;
  }

  // TODO: there's a hidden hyperparameter here
  distance += abs(target.size() - candidate.size());

  return distance;
}

int scoring::compute_wire_lengths(const connections_t &conns) {
  int lens = 0;

  // we need to compute the lengths for each individual wire i.e. unit output
  vector<vector<int>> outgoing_conns =
      compute_output_mapping_from_connections(conns);

  // we need to consider each wire individually
  for(int unit_id = 0; unit_id < UNIT_COUNT; ++ unit_id) {
    auto & out_conns = outgoing_conns[unit_id];
    if(! out_conns.empty()) {
      // make a copy and store all points of the wire, including source
      vector<int> wire(outgoing_conns[unit_id]);
      wire.push_back(unit_id);

      lens += compute_one_wire_length(wire);
    }
  }

  return lens;
}

int scoring::compute_one_wire_length(const vector<int> &wire) {
  int row_low = numeric_limits<int>::max();
  int row_high = numeric_limits<int>::lowest();

  // find longest vertical strip
  for(int unit_id : wire) {
    int unit_row = unit_id / UNIT_COLL_COUNT;

    row_low = min(row_low, unit_row);
    row_high = max(row_high, unit_row);
  }

  // store each unit in its distinct group that we unify as we go along
  utils::disj_sets groups(wire.size());

  // unify all units within a distance of 1
  for(int i = 0; i < wire.size() - 1; ++ i) {
    int row1 = wire[i] / UNIT_COLL_COUNT;
    int col1 = wire[i] % UNIT_COLL_COUNT;
    for(int j = i + 1; j < wire.size(); ++ j) {
      int row2 = wire[j] / UNIT_COLL_COUNT;
      int col2 = wire[j] % UNIT_COLL_COUNT;

      if(abs(row1 - row2) + abs(col1 - col2) == 1) {
        groups.merge(i, j);
      }
    }
  }

  // store overall wire length estimation
  int len = 0;

  // create a vertical wire as the backbone
  len += row_high - row_low;

  // assume a vertical line through each column and get the minimum spanning wire
  int min_dist = numeric_limits<int>::max();
  for(int coll_id = 0; coll_id < UNIT_COLL_COUNT; ++ coll_id) {
    int coll_dist = 0;

    // store minimum distance from a group of 1-connected units to the vertical wire
    vector<int> group_distance(wire.size(), numeric_limits<int>::max());

    for(int wid = 0; wid < wire.size(); ++ wid) {
      int unit_id = wire[wid];
      // compute distance from current unit to the vertical wire
      int unit_coll = unit_id % UNIT_COLL_COUNT;

      int group_id = groups.get_representative(wid);
      group_distance[group_id] = min(group_distance[group_id], abs(unit_coll - coll_id));
    }

    // go through all groups and accumulate minimum distances
    for(int wid = 0; wid < wire.size(); ++ wid) {
      int group_id = groups.get_representative(wid);
      if(group_distance[group_id] > 0) { // current unit group is not 1-coonected to the vertical wire
        if(group_id == wid) {
          /* This is the "leader" of the current group.
           * It has no special meaning, but we can use it to avoid connecting
           * the group twice to the vertical wire since the leader is unique
           */
          coll_dist += group_distance[group_id];
        } else {
          /* This is a unit that is 1-connected to the vertical wire.
           * Now we add its connection
           */
          coll_dist += 1;
        }
      }
    }

    // update minimum spanning distance
    min_dist = min(min_dist, coll_dist);
  }

  return len;
}
