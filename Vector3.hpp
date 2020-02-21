#ifndef VECTOR3_HPP
#define VECTOR3_HPP

#include <cmath>
#include <cassert>
#include <iostream>

#include "Constants.hpp"

template <typename T>
class Vector3
{

public:

    T x = T();
    T y = T();
    T z = T();

    Vector3(void) : x(T(0)), y(T(0)), z(T(0)) {}
    Vector3(T k) : x(k), y(k), z(k) {}
    Vector3(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}

    Vector3 & normalize(void)
    {
        T norm2 = length2();
        if (norm2 > 0) {
            T invNorm = 1. / sqrt(norm2);
            x *= invNorm;
            y *= invNorm;
            z *= invNorm;
        }

        return *this;
    }

    bool isNull(void) const
    {
        return (fabs(x) < Constants::kEpsilon &&
                fabs(y) < Constants::kEpsilon &&
                fabs(z) < Constants::kEpsilon);
    }

    T length2(void) const { return pow(x, 2) + pow(y, 2) + pow(z, 2); }
    float length(void) const { return sqrt(length2()); }

    T dot(const Vector3<T> & v) const { return Vector3<T>(x * v.x + y * v.y +  z * v.z); }

    bool hasNaN(void) const { return std::isnan(x) || std::isnan(y) || std::isnan(z); }

    // Vector3<T> operator * (const T & c) const { return Vector3<T>(x * c, y * c, z * c); }
    // Vector3<T> operator * (const Vector3<T> & v) const { return Vector3<T>(x * v.x, y * v.y, z * v.z); }

    // Vector3<T> operator = (Vector3<T> v) { return *this;}
    friend Vector3<T> operator * (const T & c, const Vector3<T> & v) { return Vector3<T>(v.x * c, v.y * c, v.z * c); }
    friend Vector3<T> operator * (const Vector3<T> & v, const T & c) { return Vector3<T>(v.x * c, v.y * c, v.z * c); }

    friend Vector3<T> operator + (const Vector3<T> & v, const T & c) { return Vector3<T>(v.x + c, v.y + c, v.z + c); }
    friend Vector3<T> operator - (const Vector3<T> & v, const T & c) { return Vector3<T>(v.x - c, v.y - c, v.z - c); }

    friend Vector3<T> operator + (const T & c, const Vector3<T> & v) { return Vector3<T>(c + v.x, c + v.y, c + v.z); }
    friend Vector3<T> operator - (const T & c, const Vector3<T> & v) { return Vector3<T>(c - v.x, c - v.y, c - v.z); }

    friend Vector3<T> operator * (const Vector3<T> & v1, const Vector3<T> & v2) { return Vector3<T>(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z); }

    friend Vector3<T> operator + (const Vector3<T> & v1, const Vector3<T> & v2) { return Vector3<T>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z); }
    friend Vector3<T> operator - (const Vector3<T> & v1, const Vector3<T> & v2) { return v1 + (-v2); }

    bool operator== (const Vector3<T> v2) const
    {
        if ((x == v2.x) && (y == v2.y) && (z == v2.z)) return true;
        return false;
    }

    bool operator!= (const Vector3<T> v2) const { return !(*this == v2); }

    Vector3<T> operator /= (const float s)
    {
        assert(s != 0);

        float inv = 1. / s;

        return Vector3<T>(x * inv, y * inv, z * inv);
    }

    Vector3<T> & operator += (const Vector3<T> & v) { x += v.x; y += v.y; z += v.z; return *this; }
    Vector3<T> & operator -= (const Vector3<T> & v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
    Vector3<T> & operator *= (const Vector3<T> & v) { x *= v.x; y *= v.y; z *= v.z; return *this; }
    Vector3<T> & operator /= (const Vector3<T> & v) { x /= v.x; y /= v.y; z /= v.z; return *this; }

    Vector3 operator - (void) const { return Vector3(-x, -y, -z); }

    friend std::ostream & operator << (std::ostream & os, const Vector3<T> & v)
    {
        os << "[" << v.x << ", " << v.y << ", " << v.z << "]";

        return os;
    }

};

#endif
