#include "../extern/catch.hpp"

#include <iostream>
#include "../src/stochastic_search.h"
#include "../src/scoring.h"

using namespace std;
using namespace scoring;

TEST_CASE("Can compute distance between polynomials", "[scoring]" ) {
  poly_t p1{3, 2, 1};
  poly_t p_empty;
  poly_t p2{3};
  poly_t p3{3, 2};
  poly_t p4{5, 4, 2};

  double d1e = compute_poly_distance(p1, p_empty);
  double d11 = compute_poly_distance(p1, p1);
  double d12 = compute_poly_distance(p1, p2);
  double d13 = compute_poly_distance(p1, p3);
  double d14 = compute_poly_distance(p1, p4);

  REQUIRE(d11 == 0.0);
  REQUIRE(d1e > d12);
  REQUIRE(d1e > d13);
  REQUIRE(d1e > d14);
  REQUIRE(d14 >= d13);
  REQUIRE(d12 > d13);
}

TEST_CASE("Can compute length for one wire #1", "[scoring]") {
  /* Implement this wire example (should be 6):
   *    012
   *    ---
   * 0: 001
   * 1: 100
   * 2: 100
   * 3: 001
   */
  vector<int> wire1 = {
    // from (row, col) to unit_id
    0 * 3 + 2,
    1 * 3 + 0,
    2 * 3 + 0,
    3 * 3 + 2,
  };

  int actual_len = compute_one_wire_length(wire1);

  REQUIRE(actual_len == 6);
}

TEST_CASE("Can compute length for one wire #2", "[scoring]") {
  /* Implement this wire example (should be 4):
   *    012
   *    ---
   * 0: 001
   * 1: 000
   * 2: 100
   */
  vector<int> wire1 = {
    // from (row, col) to unit_id
    0 * 3 + 2,
    2 * 3 + 0,
  };

  int actual_len = compute_one_wire_length(wire1);

  REQUIRE(actual_len == 4);
}

TEST_CASE("Can compute length for one wire #3", "[scoring]") {
  /* Implement this wire example (should be 4):
   *    012
   *    ---
   * 0: 001
   * 1: 010
   * 2: 000
   * 3: 001
   */
  vector<int> wire1 = {
    // from (row, col) to unit_id
    0 * 3 + 2,
    1 * 3 + 1,
    3 * 3 + 2,
  };

  int actual_len = compute_one_wire_length(wire1);

  REQUIRE(actual_len == 4);
}
