#include "stochastic_search.h"

#include <iostream>
#include <limits>
#include <algorithm>
#include <deque>

using namespace std;

// use an anonymous namespace for free functions to avoid accidental name aliasing
namespace {
void sort_canonical(vector<int> * p) {
  sort(p->begin(), p->end(), std::greater<int>());
}

/* Computes the output of a unit given its inputs that can be polynomials or invalid.
 * unit_type can be 0 (adder), 1 (multiplier) or 2 (divider)
 */
UnitOutput comput_one_unit_output(int unit_type, UnitOutput const & in1, UnitOutput const & in2) {
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


  return {true, true, poly};
}
}

StochasticSearch::StochasticSearch(const vector<int> &polynomial, int walker_count, ScoringParams params)
  : random_generator(random_device{}()),
    dist_walkers(0, walker_count - 1),
    params(params),
    poly(polynomial) {
  initialize_walkers(walker_count);

  // make sure input polynomial is in canonical form i.e. higher powers at front
  sort_canonical(& poly);
}

void StochasticSearch::train(int iteration_count, int cycle_count, int clone_count) {
  for(int iter_id = 0; iter_id < iteration_count; ++ iter_id) {
    cout << "Performing iteration [ " << iter_id + 1
         << " / " << iteration_count << " ]"
         << endl;
    for(int cycle_id = 0; cycle_id < cycle_count; ++ cycle_id) {
      perform_cycle(iter_id, cycle_id, clone_count);
    }

    // inject random noise into walkers
    inject_noise();
  }
}

void StochasticSearch::initialize_walkers(int walker_count) {
  // initialize wire connections to nil
  connections_t empty_walker(CONN_INPUT_COUNT, -1);

  // make sure we can reuse the same object by eliminating previous state
  walkers.clear();
  walkers.reserve(walker_count);
  walkers.resize(walker_count, empty_walker);
}

void StochasticSearch::perform_cycle(int iteration_id, int cycle_id, int clone_count) {
  // first compute the scores of each walker
  vector<double> scores(walkers.size(), numeric_limits<double>::lowest());

  for(int wid = 0; wid < walkers.size(); ++ wid) {
    scores[wid] = compute_score(wid);
  }

  // now perform the cloning
  int clones_performed = 0;

  while(clones_performed < clone_count) {
    int wid1 = get_random_walker_id();
    int wid2 = get_random_walker_id();

    if(wid1 != wid2) {
      if(scores[wid1] > scores[wid2]) {
        // clone walker wid1 into wid2
        walkers[wid2] = walkers[wid1];
      } else {
        // the reverse
        walkers[wid1] = walkers[wid2];
      }

      ++ clones_performed;
    }
  }
}

double StochasticSearch::compute_score(int walker_id) {
  connections_t const & walker = walkers[walker_id];

  double score = 0;

  // score having a wire connection from the input of the array
  for(int input_id = 0; input_id < walker.size(); ++ input_id) {
    if(walker[input_id] == ARRAY_INPUT_ID) {
      score += params.input_recovered_factor;
      break;
    }
  }

  // score number of units that have both inputs connected
  int count_one_input_connected = 0;
  int count_both_inputs_connected = 0;
  for(int input_id = 0; input_id < walker.size(); input_id += 2) {
    if(walker[input_id] != -1) {
      if(walker[input_id + 1] != -1) {
        ++ count_both_inputs_connected;
      } else {
        ++ count_one_input_connected;
      }
    }
  }
  if(count_both_inputs_connected > 0) {
    score += 1.0 + count_both_inputs_connected * params.unit_both_inputs_factor;
  }

  // penalize unit outputs that have a coefficient (ex: by adding x^2 and x^2 we get 2*x^2) or negative powers
  auto unit_outputs = compute_unit_outputs(walker_id);

  // score distance between unit outputs and function terms

  // score all terms that were recovered

  // extra score if the whole function is recovered (useful ?)

  // score speed prior i.e. all wire lengths

  return score;
}

unit_outputs_t StochasticSearch::compute_unit_outputs(int walker_id) {
  unit_outputs_t unit_outputs(CONN_UNIT_COUNT, {false, false, {}});

  // shorthand
  auto & walker = walkers[walker_id];

  /* Do a forward traversal starting from the array input and propagate its
   * signal to all connections. Then do the same for all units that have both
   * inputs connected.
   */

  // first construct a reverse mapping from unit_id to list of units it is connected to
  vector<vector<int>> outgoing_conns(CONN_UNIT_COUNT);

  for(int uw_id = 0; uw_id < walker.size(); ++ uw_id) {
    int in_unit_id = walker[uw_id];
    if(in_unit_id != -1) {
      int unit_id = uw_id / 2;
      outgoing_conns[in_unit_id].push_back(unit_id);
    }
  }

  // now do the propagation, starting from the input of the array because
  // it always outputs the polynomial "x"
  deque<int> propagation_front;
  propagation_front.push_back(ARRAY_INPUT_ID);
  unit_outputs[ARRAY_INPUT_ID].has_output = true;
  unit_outputs[ARRAY_INPUT_ID].is_valid = true;
  unit_outputs[ARRAY_INPUT_ID].poly = {1};

  while(! propagation_front.empty()) {
    int unit_id = propagation_front.front();
    propagation_front.pop_front();

    // compute output of current unit
    if(unit_id != ARRAY_INPUT_ID) {
      int in_unit_id1 = walker[unit_id * 2];
      int in_unit_id2 = walker[unit_id * 2 + 1];
      int unit_type = unit_id % 3;
      unit_outputs[unit_id] = comput_one_unit_output(
        unit_type,
        unit_outputs[in_unit_id1],
        unit_outputs[in_unit_id2]
      );
    }

    // propagate to its downstream units
    for(int downstream_unit_id : outgoing_conns[unit_id]) {
      // check if both inputs are connected and add it to the list to be processed
      int down_unit_in_id1 = walker[downstream_unit_id * 2];
      int down_unit_in_id2 = walker[downstream_unit_id * 2 + 1];

      // sanity check
      if(down_unit_in_id1 != unit_id && down_unit_in_id2 != unit_id) {
        cerr << "ERROR: propagation graph structure is broken for unit " << unit_id << endl;
      }

      if(down_unit_in_id1 != -1 && down_unit_in_id2 != -1) {
        // both inputs connected
        // make sure both inputs have signal flowing through
        if(unit_outputs[down_unit_in_id1].has_output &&
           unit_outputs[down_unit_in_id2].has_output) {
          propagation_front.push_back(downstream_unit_id);
        }
      }
    }

  }

  return unit_outputs;
}

void StochasticSearch::inject_noise() {
  // TODO: implement
}

int StochasticSearch::get_random_walker_id() {
  return dist_walkers(random_generator);
}
