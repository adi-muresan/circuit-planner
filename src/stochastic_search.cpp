#include "stochastic_search.h"

#include <iostream>
#include <limits>

using namespace std;

StochasticSearch::StochasticSearch(const vector<int> &polynomial, int walker_count, ScoringParams params)
  : random_generator(random_device{}()),
    dist_walkers(0, walker_count - 1),
    params(params),
    poly(polynomial) {
  initialize_walkers(walker_count);
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

  // score having a wire connection to the output of the array
  if(walker[ARRAY_OUTPUT_ID] != -1) {
    score += params.output_recovered_factor;
  }

  // score number of units that have both inputs connected
  int count_one_input_connected = 0;
  int count_both_inputs_connected = 0;
  for(int input_id = 0; input_id < ARRAY_OUTPUT_ID; input_id += 2) {
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

  // score distance between unit outputs and function terms

  // score all terms that were recovered

  // extra score if the whole function is recovered (useful ?)

  // score speed prior i.e. all wire lengths

  return score;
}

unit_outputs_t StochasticSearch::compute_unit_outputs(int walker_id) {
  // TODO: implement
  unit_outputs_t outputs(CONN_UNITS_COUNT - 1, {false, false, 0});

  return outputs;
}

void StochasticSearch::inject_noise() {
  // TODO: implement
}

int StochasticSearch::get_random_walker_id() {
  return dist_walkers(random_generator);
}
