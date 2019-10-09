#ifndef STOCHASTICSEARCH_H
#define STOCHASTICSEARCH_H

#include "definitions.h"
#include "scoring.h"
#include <vector>
#include <random>

struct NoiseParams {
  double starting_inputs_change_fraction;
  double inputs_change_decay;
  double min_inputs_change_fraction;
  double probability_change_valid_input;
  int retries_on_cycle;
};

class StochasticSearch {
  // random engine
  std::mt19937 random_generator;
  std::uniform_int_distribution<int> dist_walkers;
  std::uniform_int_distribution<int> dist_inputs;

  // hyperparameters for the scoring function i.e. tradeoffs between the
  // different score components
  scoring::ScoringParams params;

  // polynomial function to recover
  std::vector<int> poly;

  // the population of walkers
  std::vector<connections_t> walkers;

  void initialize_walkers(int walker_count);
  void perform_cycle(int iteration_id, int cycle_id, int clone_count);

  /* Implements the scoring metric that we use to drive the stochastic search.
   *
   */
  double compute_score(int walker_id);

  /* Implement graph traversal to compute what outputs each unit generates.
   */
  unit_outputs_t compute_unit_outputs(const connections_t &conns);

  // injects random noise into walkers, disallowing cycles
  void inject_noise(double iter_fraction, NoiseParams const & noise_cfg);
  void try_connect(connections_t * cconns, int input_id, std::vector<int> const & unit_ids, int retries_on_cycle);

  // utility functions for random sampling
  int get_random_walker_id();
  int get_random_input_id();

public:
  StochasticSearch(std::vector<int> const & polynomial, int walker_count, scoring::ScoringParams params);
  void train(int iteration_count, int cycle_count, int clone_count, const NoiseParams &noise);
};

#endif // STOCHASTICSEARCH_H
