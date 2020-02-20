#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <cmath>
#include <algorithm>

#include "Vector3.hpp"

namespace Utilities
{

template <typename T>
std::pair<bool, std::pair<T, T> > solveQuadratic(const T & a,
                                                 const T & b,
                                                 const T & c)
{
    T desc = b * b - 4 * a * c;
    if (desc < 0) return std::make_pair(false, std::make_pair(0, 0));
    else if (desc == 0){
        T x = -0.5 * b / a;
        return std::make_pair(true, std::make_pair(x, x));
    }

    T q = (b > 0) ?
        -0.5 * (b + sqrt(desc)):
        -0.5 * (b - sqrt(desc));

    T x1 = q / a;
    T x2 = c / q;

    // So we have it in the correct order for t0, t1
    if (x1 > x2) return std::make_pair(true, std::make_pair(x2, x1));
    else return std::make_pair(true, std::make_pair(x1, x2));
}

template <typename T>
bool solveQuadratic(const T & a,
                    const T & b,
                    const T & c,
                          T & t0,
                          T & t1)
{
    T desc = b * b - 4 * a * c;
    if (desc < 0) return false;
    else if (desc == 0){
        T x = -0.5 * b / a;
        t0 = x;
        t1 = x;
        return true;
    }

    T q = (b > 0) ?
        -0.5 * (b + sqrt(desc)):
        -0.5 * (b - sqrt(desc));

    t0 = q / a;
    t1 = c / q;

    // So we have it in the correct order for t0, t1
    if (t0 > t1) std::swap(t0, t1);
    return true;
}


template <typename T>
float dot(const Vector3<T> & v1, const Vector3<T> & v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <typename T>
float absDot(const Vector3<T> & v1, const Vector3<T> & v2)
{
    return fabs(dot(v1, v2));
}

template <typename T>
float cosAngle(const Vector3<T> & v1, const Vector3<T> & v2)
{
    return dot(v1, v2) / (v1.length * v2.length);
}

template <typename T>
float angle(const Vector3<T> & v1, const Vector3<T> & v2)
{
    return acos(cosAngle(v1, v2));
}


template <typename T>
Vector3<T> cross(const Vector3<T> & v1, const Vector3<T> & v2)
{
    float x = v1.y * v2.z - v1.z * v2.y;
    float y = v1.z * v2.x - v1.x * v2.z;
    float z = v1.x * v2.y - v1.y * v2.x;

    return Vector3<T>(x, y, z);
}

} // namespace Utilities

#endif
