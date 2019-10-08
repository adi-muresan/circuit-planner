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
  // initialize wire connections to false
  connections_t empty_walker(CONN_INPUT_SIZE, vector<bool>(CONN_OUTPUT_SIZE, false));

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

  // score having a wire connection from the input of the array

  // score having a wire connection to the output of the array

  // score number of units that have both inputs connected

  // score all terms that were recovered

  // extra score if the whole function is recovered (useful ?)

  // score speed prior i.e. all wire lengths

  return 0;
}

void StochasticSearch::inject_noise() {
  // TODO: implement
}

int StochasticSearch::get_random_walker_id() {
  return dist_walkers(random_generator);
}
