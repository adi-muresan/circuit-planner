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
UnitOutput compute_one_unit_output(int unit_type, UnitOutput const & in1, UnitOutput const & in2);

/* Compute a mapping from unit output to the units it connects to
 */
std::vector<std::vector<int>> compute_output_mapping_from_connections(connections_t const & conn);

/* Traverse the connection graph upstream and return true if there is
 * a connection from the unit with input input_id to unit_id i.e.
 * adding unit_id as a downstream connection from input_id would create a cycle.
 */
bool has_upstream_conn(const connections_t &conns, int downstream_unit_id, int upstream_unit_id);

}


#endif // PROPAGATION_H
