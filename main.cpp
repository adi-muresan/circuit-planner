#include <iostream>
#include "src/stochastic_search.h"

using namespace std;

int main(int argc, char *argv[]) {
  ScoringParams params {
    1.0,
    1.0,
    1.0
  };
  vector<int> poly {3, 7};
  StochasticSearch ss(poly, 10, params);
  ss.train(20, 30, 10);
  return 0;
}
