#include "stochastic_search.h"
#include "utils/disjoint_sets.h"
#include "propagation.h"

#include <iostream>
#include <limits>
#include <algorithm>
#include <deque>
#include <cmath>
#include <iomanip>
#include <cassert>

using namespace std;
using namespace scoring;
using namespace propagation;

StochasticSearch::StochasticSearch(const vector<int> &polynomial, int walker_count, ScoringParams params)
  : random_generator(random_device{}()),
//  : random_generator(42), // for reproducible debugging
    dist_walkers(0, walker_count - 1),
    dist_inputs(0, CONN_INPUT_COUNT - 1),
    params(params),
    poly(polynomial) {
  initialize_walkers(walker_count);

  // make sure input polynomial is in canonical form i.e. higher powers at front
  sort_canonical(& poly);
}

void StochasticSearch::train(int iteration_count, int cycle_count, int clone_count, NoiseParams const & noise_cfg) {
  for(int iter_id = 0; iter_id < iteration_count; ++ iter_id) {
    cout << "Performing iteration [ " << iter_id + 1
         << " / " << iteration_count << " ]"
         << endl;

    ScoreOutput best_score = {0, numeric_limits<double>::lowest()};

    for(int cycle_id = 0; cycle_id < cycle_count; ++ cycle_id) {
      auto score_out = perform_cycle(iter_id, cycle_id, clone_count);
      if(best_score.best_score < score_out.best_score) {
        best_score = score_out;
      }
    }

    cout << "\tbest score this iteration: "
         << setprecision(2)
         << best_score.best_score << endl
         << "\tfunction was recovered "
         << best_score.times_function_recovered
         << " times" << endl;

    // inject random noise into walkers
    double iter_fraction = (double) (iter_id + 1) / iteration_count;
    inject_noise(iter_fraction, noise_cfg);
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

ScoreOutput StochasticSearch::perform_cycle(int iteration_id, int cycle_id, int clone_count) {
  // first compute the scores of each walker
  vector<double> scores(walkers.size(), numeric_limits<double>::lowest());
  ScoreOutput best{0, numeric_limits<double>::lowest()};

  for(int wid = 0; wid < walkers.size(); ++ wid) {
    auto score_out = compute_score(wid);
    scores[wid] = score_out.best_score;

    if(best.best_score < scores[wid]) {
      best.best_score = scores[wid];
      best.times_function_recovered = score_out.times_function_recovered;
    }
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

  return best;
}

ScoreOutput StochasticSearch::compute_score(int walker_id) {
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

  // score distance between unit outputs and function terms
  auto unit_outputs = compute_unit_outputs(walker);

  vector<double> distances;
  distances.reserve(unit_outputs.size());

  for(int uid = 0; uid < unit_outputs.size(); ++ uid) {
    auto & uo = unit_outputs[uid];
    if(uo.has_output && uo.is_valid) {
      distances.push_back(compute_poly_distance(poly, uo.poly));
    }
  }

  // take the top 3 closes distances and add them to the score
  sort(distances.begin(), distances.end());

  for(int did = 0; did < 3 && did < distances.size(); ++ did) {
    // we actually want the opposite of the distance
    // take exp(-distance) because we want this to be symetrically "spikey"
    score += exp(-1.0 * distances[did]) * params.distance_factor;
  }

  // score all terms that were successfully recovered (still useful in light of the above ?)

  // extra score if the whole function is recovered by a unit output
  int times_recovered = 0;
  for(int uid = 0; uid < unit_outputs.size(); ++ uid) {
    auto & uo = unit_outputs[uid];
    if(uo.poly == poly) {
      ++ times_recovered;
    }
  }
  if(times_recovered > 0) {
    score += params.function_recovered_factor;
  }

  // score speed prior i.e. all wire lengths
  double wire_lengths = compute_wire_lengths(walker);
  score += params.speed_prior_factor * 1.0 / (1.0 + wire_lengths);

  return {times_recovered, score};
}

unit_outputs_t StochasticSearch::compute_unit_outputs(connections_t const & conns) {
  unit_outputs_t unit_outputs(CONN_UNIT_COUNT, {false, false, {}});

  /* Do a forward traversal starting from the array input and propagate its
   * signal to all connections. Then do the same for all units that have both
   * inputs connected.
   */

  // first construct a reverse mapping from unit_id to list of units it is connected to
  vector<vector<int>> outgoing_conns =
    compute_output_mapping_from_connections(conns);

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
      int in_unit_id1 = conns[unit_id * 2];
      int in_unit_id2 = conns[unit_id * 2 + 1];
      int unit_type = unit_id % 3;
      unit_outputs[unit_id] = compute_one_unit_output(
        unit_type,
        unit_outputs[in_unit_id1],
        unit_outputs[in_unit_id2]
      );
    }

    // propagate to its downstream units
    for(int downstream_unit_id : outgoing_conns[unit_id]) {
      // check if both inputs are connected and add it to the list to be processed
      int down_unit_in_id1 = conns[downstream_unit_id * 2];
      int down_unit_in_id2 = conns[downstream_unit_id * 2 + 1];

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

/* Inject some noise into all the random walkers.
 * In general it's good to inject more noise in the beginning and less towards
 * the end of training when we already have partial solutions.
 *
 * The noise is in the form of changing the input connections of a unit:
 * - give higher chances to change an input that is not connected
 * - only connect to units that have a valid output
 */
void StochasticSearch::inject_noise(double iter_fraction, const NoiseParams &noise_cfg) {
  double fraction_to_change =
      noise_cfg.starting_inputs_change_fraction *

      // exponential decay the number of inputs we change
      pow(1.0 - noise_cfg.inputs_change_decay, iter_fraction * 10);

  // cap the minimum to make sure we always change at least a few inputs
  fraction_to_change = max(fraction_to_change, noise_cfg.min_inputs_change_fraction);
  int inputs_to_change = CONN_INPUT_COUNT * fraction_to_change;

  bernoulli_distribution change_valid_input(noise_cfg.probability_change_valid_input);

  for(auto & walker : walkers) {
    auto unit_outputs = compute_unit_outputs(walker);

    // compute units with valid outputs
    vector<int> units_with_valid_outputs;
    for(int unit_id = 0; unit_id < unit_outputs.size(); ++ unit_id) {
      if(unit_outputs[unit_id].is_valid) {
        units_with_valid_outputs.push_back(unit_id);
      }
    }

    for(int cid = 0; cid < inputs_to_change; ++ cid) {
      int input_id = get_random_input_id();
      int upstream_unit_id = walker[input_id];

      if(upstream_unit_id >= 0 && unit_outputs[upstream_unit_id].is_valid) {
        // input is connected to a wire producing a valid signal
        // sample the config Bernoulli distribution to see if we should change it
        if(change_valid_input(random_generator)) {
          try_connect(& walker, input_id, units_with_valid_outputs, noise_cfg.retries_on_cycle);
        }
      } else {
        try_connect(& walker, input_id, units_with_valid_outputs, noise_cfg.retries_on_cycle);
        // TODO: could look into updating the units_with_valid_outputs on the fly here
      }
    }
  }
}

void StochasticSearch::try_connect(connections_t * conns, int input_id, const std::vector<int> &unit_ids, int retries_on_cycle) {
  assert(unit_ids.size() > 0);

  uniform_int_distribution<int> dist_units(0, unit_ids.size() - 1);
  int tries = 0;
  bool have_connected = false;
  int unit_id = input_id / 2;

  do {
    int sampled_index = dist_units(random_generator);
    int target_unit_id = unit_ids[sampled_index];

    // connect only if this would not introduce a cycle
    if(! has_upstream_conn(*conns, target_unit_id, unit_id)) {
      conns->at(input_id) = target_unit_id;
      have_connected = true;
    }

    ++ tries;
  } while(tries < retries_on_cycle && ! have_connected);
}

int StochasticSearch::get_random_walker_id() {
  return dist_walkers(random_generator);
}

int StochasticSearch::get_random_input_id() {
  return dist_inputs(random_generator);
}
