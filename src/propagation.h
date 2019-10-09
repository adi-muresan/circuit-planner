#ifndef PROPAGATION_H
#define PROPAGATION_H

#include <vector>
#include "definitions.h"

namespace propagation {
/* Sort a polynomial's powers in canonycal order */
void sort_canonical(std::vector<int> * p);

/* Computes the output of a unit given its inputs that can be polynomials or invalid.
 * unit_type can be 0 (adder), 1 (multiplier) or 2 (divider)
 */
UnitOutput comput_one_unit_output(int unit_type, UnitOutput const & in1, UnitOutput const & in2);

/* Compute a mapping from unit output to the units it connects to
 */
std::vector<std::vector<int>> compute_output_mapping_from_connections(connections_t const & conn);

}


#endif // PROPAGATION_H
