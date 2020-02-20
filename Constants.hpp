#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <limits>

namespace Constants
{

const float kEpsilon  = std::numeric_limits<float>::epsilon();
const float kInfinity = std::numeric_limits<float>::max();

const float k1_4 = 1 / 4.;

} // namespace Constants

#endif
