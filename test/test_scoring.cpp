#include "../extern/catch.hpp"

#include <iostream>
#include "../src/stochastic_search.h"

using namespace std;

TEST_CASE("Can run stochastic search", "[stochastic_search]" ) {
  // TODO: test all free functions of polynomial propagation and computation
  ScoringParams params {
    1.0,
    1.0,
    1.0,
    0.2,
    1.0,
    100.0,
    10.0,
    10.0
  };
  vector<int> poly {3, 7};
  StochasticSearch ss(poly, 10, params);
  ss.train(20, 30, 10);
}

// TODO: more tests
