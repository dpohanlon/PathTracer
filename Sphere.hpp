#ifndef SPHERE_HPP
#define SPHERE_HPP

#include <cmath>

#include "Vector3.hpp"
#include "Utilities.hpp"
#include "Ray.h"

template <typename T>
class Sphere
{
public:

    Vector3<T> center;
    T radius;
    T radius2;

    Vector3<float> surfaceColor;
    Vector3<float> emissionColor;

    T transparency;
    T reflection;

    bool emissive;

    Sphere(const Vector3<T> & c,
           const T & r,
           const Vector3<float> & sc,
           const T & trs = 0,
           const T & rfl = 0,
           const Vector3<float> & ec = 0
          ) :
            center(c),
            radius(r),
            radius2(r * r),
            surfaceColor(sc),
            emissionColor(ec),
            transparency(trs),
            reflection(rfl),
            emissive(!ec.isNull())
    {}

    // bool intersect(const Vector3<T> & rayOri,
    //                const Vector3<T> & rayDir,
    //                T * t0 = nullptr,
    //                T * t1 = nullptr) const
    // {
    //     Vector3<T> l = rayOri - center;
    //     T tca = Utilities::dot(l, rayDir);
    //
    //     if (tca < 0) return false;
    //
    //     T d2 = Utilities::dot(l, l) - tca * tca;
    //
    //     if (d2 > radius2) return false;
    //
    //     T thc = sqrt(radius2 - d2);
    //
    //     if (t0 != nullptr && t1 != nullptr) {
    //         *t0 = tca - thc;
    //         *t1 = tca + thc;
    //     }
    //
    //     return true;
    // }

    bool intersect(const Vector3<T> & rayOri,
                   const Vector3<T> & rayDir,
                   T & t0,
                   T & t1) const
    {
        Vector3<T> l = rayOri - center;

        float a = Utilities::dot(rayDir, rayDir);
        float b = 2 * Utilities::dot(rayDir, l);
        float c = Utilities::dot(l, l) - radius2;

        if (!Utilities::solveQuadratic(a, b, c, t0, t1)) {
            return false;
        } else {
            return true;
        }

    }

    // bool intersect(Ray & r) const
    // {
    //
    //     Vector3<float> l = r.origin - center;
    //     float a = Utilities::dot(r.direction, r.direction);
    //     float b = 2 * Utilities::dot(r.direction, l);
    //     float c = Utilities::dot(l, l) - radius2;
    //
    //     float t0 = 0;
    //     float t1 = 0;
    //
    //     if (!Utilities::solveQuadratic(a, b, c, t0, t1))  return false;
    //
    //     if (t0 > r.maxT || t1 < r.minT) return false;
    //     float tHit = t0;
    //     if (t0 < r.minT) {
    //         tHit = t1;
    //         if (tHit > r.maxT) return false;
    //     }
    //
    //     r.t = t0;
    //
    //     return true;
    // }

    Vector3<float> hitNormal(Vector3<float> & pHit)
    {
        Vector3<float> norm = pHit - this->center;
        norm.normalize();

        return norm;
    }

    void illuminate(Vector3<float> & pHit,
                    Vector3<float> & lightDirection,
                    Vector3<float> & lightIntensity,
                    float & distance)
    {
        lightDirection = pHit - this->center;
        distance = lightDirection.length();
        lightDirection.normalize();

        lightIntensity = this->emissionColor; // Magnitude of ec is intensity
    }

    bool is(const Sphere<T> & s) const
    {
        if (s.center != center) return false;
        if (s.radius != radius) return false;
        if (s.surfaceColor != surfaceColor) return false;
        if (s.emissionColor != emissionColor) return false;
        if (s.transparency != transparency) return false;
        if (s.reflection != reflection) return false;

        return true;
    }

};

#endif
