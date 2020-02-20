#include "Ray.h"

Ray::Ray(Point p, Vec v, float t_) : origin(p), direction(v), t(t_)
{
    // origin.normalize();
    direction.normalize();
}

Ray::Ray(Point p, Vec v) : origin(p), direction(v), t(Constants::kEpsilon)
{
    // origin.normalize();
    direction.normalize();
}

Ray::Ray(const Point & o, const Vec & d, const Ray & parent, float start, float end)
         : origin(o), direction(d), minT(start), maxT(end),
           depth(parent.depth + 1), time(parent.time)
{

}
