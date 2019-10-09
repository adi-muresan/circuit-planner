#include "propagation.h"

#include <algorithm>
#include <iostream>

using namespace std;

void propagation::sort_canonical(std::vector<int> *p) {
  sort(p->begin(), p->end(), std::greater<int>());
}

UnitOutput propagation::comput_one_unit_output(int unit_type, const UnitOutput &in1, const UnitOutput &in2) {
  if(! in1.has_output || ! in2.has_output) {
    // if an input does not have a signal flowing through it then do not propagate
    // this should never happen btw
    return {false, false, {}};
  }
  if(! in1.is_valid || ! in2.is_valid) {
    // if one of the inputs is an invalid polynomial then do not propagate it
    return {true, false, {}};
  }

  vector<int> poly;

  // implement polynomial addition, multiplication and division
  switch(unit_type) {

  case 0: // addition -- just concat all terms
    // copy of the first polynomial
    poly = in1.poly;
    for(int p : in2.poly) {
      poly.push_back(p);
    }
    break;

  case 1: // multiplication
    poly.reserve(in1.poly.size() + in2.poly.size());
    for(int p1 : in1.poly) {
      for(int p2 : in2.poly) {
        poly.push_back(p1 * p2);
      }
    }
    break;

  case 2: // division. For now only by polynomials with a single term
    if(in2.poly.size() > 1 || in2.poly.size() == 0) {
      return {true, false, {}};
    } else {
      int divider = in2.poly[0];

      // copy of the first polynomial
      poly = in1.poly;

      for(int& p : poly) {
        p -= divider;
      }
    }
    break;

  default: // this should never happen
    cerr << "ERROR: Received unit of invalid type: " << unit_type << endl;
    return {true, false, {}};
  }

  sort_canonical(& poly);

  // check for duplicate powers i.e. a polynomial of the form "2*x" which is invalid
  for(int pid = 1; pid < poly.size(); ++ pid) {
    if(poly[pid - 1] == poly[pid]) {
      return {true, false, {}};
    }
  }

  // check for negative powers
  if(poly.back() <= 0) {
    return {true, false, {}};
  }

  return {true, true, poly};
}

std::vector<std::vector<int> > propagation::compute_output_mapping_from_connections(const connections_t &conn) {
  vector<vector<int>> outgoing_conns(CONN_UNIT_COUNT);

  for(int uw_id = 0; uw_id < conn.size(); ++ uw_id) {
    int in_unit_id = conn[uw_id];
    if(in_unit_id != -1) {
      int unit_id = uw_id / 2;
      outgoing_conns[in_unit_id].push_back(unit_id);
    }
  }

  return outgoing_conns;
}

bool propagation::had_upstream_conn(const connections_t &conns, int downstream_unit_id, int upstream_unit_id) {
  // handle special case; the array input unit has no upstream units
  if(downstream_unit_id == ARRAY_INPUT_ID) {
    return false;
  }

  vector<int> queue;
  queue.push_back(downstream_unit_id);

  // standard graph traversal
  while(!queue.empty()) {
    int current_unit_id = queue.back();
    queue.pop_back();

    int upstream_unit_id1 = conns[current_unit_id * 2];
    if(upstream_unit_id1 != -1) {
      if(upstream_unit_id1 == upstream_unit_id) {
        return true;
      }
      queue.push_back(upstream_unit_id1);
    }

    int upstream_unit_id2 = conns[current_unit_id * 2 + 1];
    if(upstream_unit_id2 != -1) {
      if(upstream_unit_id2 == upstream_unit_id) {
        return true;
      }
      queue.push_back(upstream_unit_id2);
    }
  }

  return false;
}
