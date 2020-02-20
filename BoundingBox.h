#ifndef BOUNDINGBOX_H
#define BOUNDINGBOX_H

#include "Point.h"
#include "Constants.hpp"

class BoundingBox
{
public:
    BoundingBox(void);
    BoundingBox(const Point & p);
    BoundingBox(const Point & p1, const Point & p2);

    BoundingBox Union(const BoundingBox & b, const Point & p);

    ~BoundingBox();

    Point pMax;
    Point pMin;
};

#endif
