#include "../extern/catch.hpp"

#include <iostream>
#include "../src/propagation.h"
#include "../src/definitions.h"

using namespace std;
using namespace propagation;

TEST_CASE("Can find upstream connection", "[propagation]" ) {
  connections_t conns(CONN_INPUT_COUNT, -1);

  // simulate a few connections
  int unit1_id = 3;
  int unit2_id = 17;
  conns[unit1_id * 2] = ARRAY_INPUT_ID;
  conns[unit1_id * 2 + 1] = 123;

  conns[unit2_id * 2] = unit1_id;

  REQUIRE(has_upstream_conn(conns, unit2_id, unit1_id));
  REQUIRE(!has_upstream_conn(conns, unit1_id, unit2_id));

  REQUIRE(has_upstream_conn(conns, unit1_id, ARRAY_INPUT_ID));
  REQUIRE(has_upstream_conn(conns, unit2_id, ARRAY_INPUT_ID));
}

TEST_CASE("Can compute polynomial operations", "[propagation]" ) {
  UnitOutput p1{true, true, {3, 2}};
  UnitOutput p2{true, true, {5, 7}};
  vector<int> expected_add{7, 5, 3, 2};
  vector<int> expected_mult{10, 9, 8, 7};

  int add = 0;
  int multiply = 1;
  int divide = 2;

  auto output = compute_one_unit_output(add, p1, p2);
  REQUIRE(output.has_output);
  REQUIRE(output.is_valid);
  REQUIRE(output.poly == expected_add);

  output = compute_one_unit_output(multiply, p1, p2);
  REQUIRE(output.has_output);
  REQUIRE(output.is_valid);
  REQUIRE(output.poly == expected_mult);

  output = compute_one_unit_output(divide, p1, p2);
  REQUIRE(output.has_output);
  REQUIRE(! output.is_valid);

  // test poly division
  p1.poly = {3};
  vector<int> expected_div{4, 2};

  output = compute_one_unit_output(divide, p2, p1);
  REQUIRE(output.has_output);
  REQUIRE(output.is_valid);
  REQUIRE(output.poly == expected_div);
}
