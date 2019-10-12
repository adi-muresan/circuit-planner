#ifndef SCORING_H
#define SCORING_H

#include <vector>
#include "definitions.h"

namespace scoring {
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

/* Estimate a "distance" between a target polynomial and a candidate.
 * This distance is not symmetrical since our goal is to recover the target.
 * The distance is always positive or zero.
 *
 * A distance close to 0 means the candidate is close to the target.
 *
 * One way to compute the distance is for each term in the target find the
 * closest term in the candidate and use their distance.
 *
 * We also must account for the difference in number of terms in the two polynomials.
 */
double compute_poly_distance(std::vector<int> const & target, std::vector<int> const & candidate);

/* Compute lengths of all the wires, accounting for wire reuse.
 *
 * The physical structure of the array is simple and allows a minium spanning
 * wire for each type of connection. Computing the minimum spanning length is
 * not straightforward, which is why we will use a heuristic to approximate it.
 * Since the array is a lot longer along the Y direction, there will be longer
 * lines forming in this direction.
 *
 * Heuristic:
 * - add wire for any two points within a Manhattan distance of 1
 * - find the longest "vertical" span
 * - for each column construct a vertical line between with the span above
 * - connect all other points and point groups to the constructed line
 * - repeat for every column in Y, keeping the minimum length
 *
 * Example:
 *
 *    012
 *    ---
 * 0: 001
 * 1: 100
 * 2: 100
 * 3: 001
 *
 * h: 202 (histogram along Y)
 *
 * step 1: unite the two 1-neighbors of the first column
 * step 2: longest vertical span is from row 0 to row 3 i.e. a length of 3
 *
 * step 3-1: assume a vertical wire with the above span on the first column
 * step 3-2: connect remaining point groups in the third column to this line
 * step 3-3: resulting wire will have a length of 7
 *
 * step 4-1: assume a vertical wire with the above span on the third column
 * step 4-2: connect remaining point groups in the first column to this line
 * step 4-3: resulting wire will have a length of 6
 *
 * step 4: return a minimum length of 6.
 *
 *
 * Note: the current heuristic overestimates wire lengths i.e. there is no guarantee
 * that the solution here will be minimal wrt. spanning wire length.
 *
 * TODO: Find a better way of computing minimum wire lengths.
 *
 * TODO: Account for wire length from the input of the array to the first
 * unit and from the last unit to the output of the array.
 */
int compute_wire_lengths(connections_t const & conns);

// Implement logic described above for one wire
int compute_one_wire_length(std::vector<int> const & wire);

}


#endif // SCORING_H
