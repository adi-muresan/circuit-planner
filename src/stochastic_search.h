#ifndef STOCHASTICSEARCH_H
#define STOCHASTICSEARCH_H

#include <vector>
#include <random>

/* This represents the search space of possible wires and their connections.
 * Store wire connections as a 2D matrix of size 151 x 301:
 *
 * - 3 x 50 + 1 outputs representing the signal source of a wire
 * i.e. 3 unit types, one output each
 *
 * - 3 x 50 x 2 + 1 inputs representing possible wire outgoing connections
 * i.e. 3 unit types, 2 inputs each and one output for the whole array
 *
 *
 * Arithmetic units (Adder, Multiplier, Divider) are physically laid out in the
 * following way:
 *     012
 *     ---
 *  0: AMD
 *  1: AMD
 *     ...
 * 49: AMD
 *
 * connection_t[i] - the wire connecting the output of unit at position
 * (i/3, i % 3) in the physical array to the inputs of any units
 *
 * connection_t[i][j] - here "j" is the index to one of the inputs (input j % 2)
 * of the unit at position ((j / 2) / 3, (j / 2) % 3) in the physical array
 *
 * Special wires:
 * - connection_t[150] - the wire connecting the input of the array
 * - connection_t[i][300] - the connection to the output of the array
 */
typedef std::vector<std::vector<bool>> connections_t;

const std::size_t CONN_INPUT_SIZE = 3 * 50 + 1;
const std::size_t CONN_OUTPUT_SIZE = 3 * 50 * 2 + 1;

struct ScoringParams {
  double unit_single_input_penalty;
  double term_recovered_factor;
  double function_recovered_factor;
  double speed_prior_factor;
};

class StochasticSearch {
  // random engine
  std::mt19937 random_generator;
  std::uniform_int_distribution<int> dist_walkers;

  // hyperparameters for the scoring function i.e. tradeoffs between the
  // different score components
  ScoringParams params;

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

  // injects random noise into walkers, disallowing cycles
  void inject_noise();

  // utility functions for random sampling
  int get_random_walker_id();

public:
  StochasticSearch(std::vector<int> const & polynomial, int walker_count, ScoringParams params);
  void train(int iteration_count, int cycle_count, int clone_count);
};

#endif // STOCHASTICSEARCH_H
