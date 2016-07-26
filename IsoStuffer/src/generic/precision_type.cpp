#include "precision_type.hpp"
#include <limits>

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

    template<>
    const float PrecisionType<float>::EPS = (float)1E-8;

    template<>
    const float PrecisionType<float>::MA_EPS = std::numeric_limits<float>::epsilon();

    template <>
    const double PrecisionType<double>::EPS = 1E-12;
    template<>
    const double PrecisionType<double>::MA_EPS = std::numeric_limits<double>::epsilon();

#ifdef USE_NAMESPACE
}
#endif
