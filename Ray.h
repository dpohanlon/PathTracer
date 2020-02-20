#ifndef RAY_H
#define RAY_H

#include "Vec.h"
#include "Point.h"
#include "Constants.hpp"

class Ray
{
public:
    Ray(void);
    Ray(Point p, Vec v, float t);
    Ray(Point p, Vec v);
    Ray(const Point & o, const Vec & d, const Ray & parent,
        float start, float end = Constants::kInfinity);
    ~Ray() {}

    Point operator()(float t) const { return origin + t * direction; }

// private:
    Point origin;
    Vec direction;

    float t = Constants::kEpsilon;

    float minT = Constants::kEpsilon;
    float maxT = Constants::kInfinity;

    int depth = 0;
    float time = 0;

};

#endif
