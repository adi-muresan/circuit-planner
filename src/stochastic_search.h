#ifndef STOCHASTICSEARCH_H
#define STOCHASTICSEARCH_H

#include "definitions.h"
#include "scoring.h"
#include <vector>
#include <random>

class StochasticSearch {
  // random engine
  std::mt19937 random_generator;
  std::uniform_int_distribution<int> dist_walkers;

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
  unit_outputs_t compute_unit_outputs(int walker_id);

  // injects random noise into walkers, disallowing cycles
  void inject_noise();

  // utility functions for random sampling
  int get_random_walker_id();

public:
  StochasticSearch(std::vector<int> const & polynomial, int walker_count, scoring::ScoringParams params);
  void train(int iteration_count, int cycle_count, int clone_count);
};

#endif // STOCHASTICSEARCH_H
