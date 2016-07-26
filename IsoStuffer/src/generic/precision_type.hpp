#ifndef CARBINE_PRECISION_TYPE_HPP
#   define CARBINE_PRECISION_TYPE_HPP

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

    template <typename T>
    struct PrecisionType
    {
        static const T EPS;
        static const T MA_EPS;
    };

#ifdef USE_NAMESPACE
}
#endif

#endif
