#include "BoundingBox.h"

BoundingBox::BoundingBox(void)
    : pMin(Point(Constants::kInfinity, Constants::kInfinity, Constants::kInfinity)),
      pMax(Point(-Constants::kInfinity, -Constants::kInfinity, -Constants::kInfinity))
{

}

BoundingBox::BoundingBox(const Point & p) : pMin(p), pMax(p)
{

}

BoundingBox::BoundingBox(const Point & p1, const Point & p2)
{
    pMin = Point(fmin(p1.x, p2.x), fmin(p1.y, p2.y), fmin(p1.z, p2.z));
    pMax = Point(fmax(p1.x, p2.x), fmax(p1.y, p2.y), fmax(p1.z, p2.z));
}