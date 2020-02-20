#ifndef POINT_H
#define POINT_H

#include "Vector3.hpp"

class Point : public Vector3<float>
{
public:

    Point(void) : Vector3<float>() {}
    Point(float c) : Vector3(c) {}
    Point(float x, float y, float z) : Vector3(x, y, z) {}
    Point(Vector3<float> v) : Vector3(v) {}

    float weight = 1.0;

};

#endif
