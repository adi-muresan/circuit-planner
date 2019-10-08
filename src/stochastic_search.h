#ifndef STOCHASTICSEARCH_H
#define STOCHASTICSEARCH_H

#include <vector>
#include <random>

/* This represents the search space of possible wires and their connections.
 * Store wire connections as a 2D matrix of size 151 x 301:
 *
 * - inputs for the same unit are consecutive and start at a multiple of 2
 * - we do not have a wire for the connection to the output of the array
 * because we'll just assume a connection is made when one of the arithmetic
 * units manages to recover the desired function
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

// TODO: revise the above
typedef std::vector<int> connections_t;

const std::size_t UNIT_COUNT = 3 * 50;
const std::size_t CONN_UNIT_COUNT = UNIT_COUNT + 1;
const std::size_t CONN_INPUT_COUNT = 3 * 50 * 2;

// last element represents the input of the physical array
const std::size_t ARRAY_INPUT_ID = UNIT_COUNT - 1;

/* A unit can generate an output, but it does not have to be always valid.
 * Ex: x^2 + x^2 = 2*x^2, which we do not want to propagate further
 * Ex: x / (x + 1) is not a valid polynomial
 */
struct UnitOutput {
  // true if the unit produces an output i.e. has inputs that have a signal flowing through
  bool has_output;

  // true if the output is valid i.e. is a polynomial of the type we expect
  bool is_valid;

  // the actual polynomial the current unit is outputting
  std::vector<int> poly;
};

typedef std::vector<UnitOutput> unit_outputs_t;

struct ScoringParams {
  double input_recovered_factor;
  double output_recovered_factor;
  double unit_single_input_penalty;
  double unit_both_inputs_factor;
  double term_recovered_factor;
  double function_recovered_factor;
  double distance_factor;
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

  /* Implement graph traversal to compute what outputs each unit generates.
   */
  unit_outputs_t compute_unit_outputs(int walker_id);

  // injects random noise into walkers, disallowing cycles
  void inject_noise();

  // utility functions for random sampling
  int get_random_walker_id();

public:
  StochasticSearch(std::vector<int> const & polynomial, int walker_count, ScoringParams params);
  void train(int iteration_count, int cycle_count, int clone_count);
};

#endif // STOCHASTICSEARCH_H
