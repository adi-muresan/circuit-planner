#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <vector>

/* This represents the search space of possible wires and their connections.
 * Store wire connections as an array of size 301, one entry per unit input:
 *
 * - 3 x 50 x 2 inputs representing possible wire outgoing connections
 * i.e. 3 unit types, 2 inputs each
 *
 * - inputs for the same unit are consecutive and start at a multiple of 2
 *
 * - we do not have a wire for the connection to the output of the array
 * because we'll just assume a connection is made when one of the arithmetic
 * units manages to recover the desired function
 *
 * - the value stored in the array is the ID of the unit connected to its input,
 * there are 3 x 50 units + 1 for the input of the array
 *
 * - the ID 150 is special and represents tha output of the input unit i.e.
 * the input of the whole array, hardwired to the polynomial "x"
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
 * connection_t[i] - the unit connected to the input "i % 2" of unit at position
 * ("(i / 2) / 3", "(i / 2) % 3") in the physical array to the inputs of any units
 */
typedef std::vector<int> connections_t;

const std::size_t UNIT_ROW_COUNT = 50;
const std::size_t UNIT_COLL_COUNT = 3;
const std::size_t UNIT_COUNT = UNIT_ROW_COUNT * UNIT_COLL_COUNT;
const std::size_t CONN_UNIT_COUNT = UNIT_COUNT + 1;
const std::size_t CONN_INPUT_COUNT = 3 * 50 * 2;

// last element represents the input of the physical array
const std::size_t ARRAY_INPUT_ID = CONN_UNIT_COUNT - 1;

typedef std::vector<int> poly_t;

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

#endif // DEFINITIONS_H
